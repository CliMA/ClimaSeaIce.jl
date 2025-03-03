using Oceananigans.Utils: Time
using Oceananigans.Fields: flattened_unique_values
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!

using ClimaSeaIce.SeaIceMomentumEquations: step_momentum!

import Oceananigans.Models: update_model_field_time_series!

const FESeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:ForwardEulerTimeStepper}

function time_step!(model::FESeaIceModel, Δt; callbacks = [])
    
    # Be paranoid and update state at iteration 0
    model.clock.iteration == 0 && update_state!(model)

    compute_tendencies!(model, Δt)
    step_tracers!(model, Δt)

    # TODO: This is an implicit (or split-explicit) step to advance momentum!
    step_momentum!(model, model.dynamics, Δt, 1)

    tick!(model.clock, Δt)
    update_state!(model)

    return nothing
end

function step_tracers!(model::SIM, Δt)
    grid = model.grid
    arch = architecture(grid)

    h    = model.ice_thickness
    ℵ    = model.ice_concentration
    hmin = model.ice_consolidation_thickness
    tracers = model.tracers

    Gⁿ = model.timestepper.Gⁿ
    
    launch!(arch, grid, :xy, _step_tracers!, h, ℵ, hmin, tracers, Gⁿ, Δt)

    return nothing
end

# Thickness and concentration are updated
# We compute hⁿ⁺¹ and ℵⁿ⁺¹ in the same kernel to account for ridging: 
# if ℵ > 1, we reset the concentration to 1 and adjust the thickness 
# to conserve the total ice volume in the cell.
@kernel function _step_tracers!(h, ℵ, hmin, tracers, Gⁿ, Δt)
    i, j = @index(Global, NTuple)
    k = 1
    
    Ghⁿ = Gⁿ.h
    Gℵⁿ = Gⁿ.ℵ
    
    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁻ = hmin[i, j, k]
        h⁺ = h[i, j, k] + Δt * Ghⁿ[i, j, k]
        ℵ⁺ = ℵ[i, j, k] + Δt * Gℵⁿ[i, j, k]

        ℵ⁺ = max(zero(ℵ⁺), ℵ⁺) # Concentration cannot be negative, clip it up
        h⁺ = max(zero(h⁺), h⁺) # Thickness cannot be negative, clip it up

        ht, ℵt = conservative_adjustment(h⁺, h⁻, ℵ⁺)

        ℵ[i, j, k] = ℵt
        h[i, j, k] = ht
    end 
end

# If h < hmin we reset the thickness to h⁻ and adjust the concentration accordingly
# to maintain a constant ice volume. 
# A no ice condition is represented by h = hmin and ℵ = 0 since ice_volume = (h * ℵ)
# The thickness should _NEVER_ be zero! 
@inline function conservative_adjustment(h⁺, h⁻, ℵ⁺)

    # Remove ice if h⁺ == 0
    thin_ice = (0 < h⁺ < h⁻) # Thin ice condition

    ht = ifelse(thin_ice, h⁻, h⁺) # Non conservative adjustement of thickness
    V⁺ = h⁺ * ℵ⁺ # Total sea ice volume
    ht = ifelse(ℵ⁺ > 1, ht * ℵ⁺, ht)
    ℵt = ifelse(ht == 0, zero(ℵ⁺), ℵ⁺)
    ht = ifelse(ht == 0, h⁻, ht)
    ℵt = ifelse(ℵt > 1, one(ℵt), ℵt)

    return ht, ℵt
end

function update_state!(model::SIM)
    
    foreach(prognostic_fields(model)) do field
        mask_immersed_field!(field)
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
