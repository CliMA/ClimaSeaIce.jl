using Oceananigans.Utils: Time
using Oceananigans.Fields: flattened_unique_values
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!
using Oceananigans.ImmersedBoundaries: mask_immersed_field_xy!

using ClimaSeaIce.SeaIceMomentumEquations: time_step_momentum!
using ClimaSeaIce.SeaIceThermodynamics: thermodynamic_time_step!

import Oceananigans.Models: update_model_field_time_series!

const FESeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:ForwardEulerTimeStepper}

# We separate the thermodynamic step from the advection (dynamic) step.
# The thermodynamic step is column physics and is performed all at once.
function time_step!(model::FESeaIceModel, Δt; callbacks = [])
    
    # Be paranoid and update state at iteration 0
    model.clock.iteration == 0 && update_state!(model)

    # Perform the thermodynamic step
    thermodynamic_time_step!(model, model.ice_thermodynamics, Δt)

    # Compute advective tendencies and update 
    # advected tracers
    compute_tendencies!(model, Δt)
    dynamic_time_step!(model, Δt)

    # This is an implicit (or split-explicit) step to advance momentum.
    time_step_momentum!(model, model.dynamics, Δt)

    tick!(model.clock, Δt)
    update_state!(model)

    return nothing
end

function dynamic_time_step!(model::SIM, Δt)
    grid = model.grid
    arch = architecture(grid)

    h = model.ice_thickness
    ℵ = model.ice_concentration
    tracers = model.tracers

    Gⁿ = model.timestepper.Gⁿ
    
    launch!(arch, grid, :xy, _dynamic_step_tracers!, h, ℵ, tracers, Gⁿ, Δt)

    return nothing
end

# Thickness and concentration are updated
# We compute hⁿ⁺¹ and ℵⁿ⁺¹ in the same kernel to account for ridging: 
# if ℵ > 1, we reset the concentration to 1 and adjust the thickness 
# to conserve the total ice volume in the cell.
@kernel function _dynamic_step_tracers!(h, ℵ, tracers, Gⁿ, Δt)
    i, j = @index(Global, NTuple)
    k = 1
    
    Ghⁿ = Gⁿ.h
    Gℵⁿ = Gⁿ.ℵ

    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁺ = h[i, j, k] + Δt * Ghⁿ[i, j, k]
        ℵ⁺ = ℵ[i, j, k] + Δt * Gℵⁿ[i, j, k]

        ℵ⁺ = max(zero(ℵ⁺), ℵ⁺) # Concentration cannot be negative, clip it up
        h⁺ = max(zero(h⁺), h⁺) # Thickness cannot be negative, clip it up

        # Ridging and rafting caused by the advection step
        V⁺ = h⁺ * ℵ⁺
        
        ℵ[i, j, k] = ifelse(ℵ⁺ > 1, one(ℵ⁺), ℵ⁺)
        h[i, j, k] = ifelse(ℵ⁺ > 1, V⁺, h⁺)

        advance_tracers!(tracers, i, j, k, Gⁿ, Δt)
    end 
end

advance_tracers!(::EmptyTuples, args...) = nothing

function advance_tracers!(tracers, i, j, k, G, Δt)
    # Assumption! The tracer tendencies are the first ones
    for n in eachindex(tracers)
        _advance_tracer!(tracers[n], i, j, k, G[n], Δt)
    end
end

_advance_tracer!(tracer, i, j, k, G, Δt) = @inbounds tracer[i, j, 1] += Δt * G[i, j, k]
_advance_tracer!(::ConstantField, i, j, k, G, Δt) = nothing
_advance_tracer!(::ZeroField, i, j, k, G, Δt) = nothing
_advance_tracer!(::OneField, i, j, k, G, Δt) = nothing
_advance_tracer!(::Nothing, i, j, k, G, Δt) = nothing

function update_state!(model::SIM)
    
    foreach(prognostic_fields(model)) do field
        mask_immersed_field_xy!(field, k=size(model.grid, 3))
    end

    fill_halo_regions!(prognostic_fields(model), model.clock, fields(model))

    update_model_field_time_series!(model, model.clock)

    return nothing
end

function update_model_field_time_series!(model::SeaIceModel, clock::Clock)
    time = Time(clock.time)

    possible_fts = (model.tracers, model.external_heat_fluxes, model.dynamics)
    time_series_tuple = extract_field_time_series(possible_fts)
    time_series_tuple = flattened_unique_values(time_series_tuple)

    for fts in time_series_tuple
        update_field_time_series!(fts, time)
    end

    return nothing
end
