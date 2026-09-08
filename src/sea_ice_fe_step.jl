using Oceananigans.Fields: flattened_unique_values, ZeroField
using Oceananigans.ImmersedBoundaries: mask_immersed_field_xy!
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!
using Oceananigans.Units: Time

using .SeaIceDynamics: time_step_momentum!
using .SeaIceThermodynamics: thermodynamic_time_step!

const FESeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:Any, <:ForwardEulerTimeStepper}

# We separate the thermodynamic step from the advection (dynamic) step.
# The thermodynamic step is column physics and is performed all at once.
function Oceananigans.TimeSteppers.time_step!(model::FESeaIceModel, Δt; kwargs...)

    # Be paranoid and update state at iteration 0
    model.clock.iteration == 0 && update_state!(model)

    # Compute advective tendencies and update advected tracers
    compute_tendencies!(model, Δt)

    # This is an implicit (or split-explicit) step to advance momentum.
    time_step_momentum!(model, model.dynamics, Δt)

    # Dynamic step for tracers
    dynamic_time_step!(model, Δt)

    # Perform the thermodynamic step
    thermodynamic_time_step!(model, model.ice_thermodynamics, model.snow_thermodynamics, Δt)

    tick!(model.clock, Δt)
    update_state!(model)

    return nothing
end

dynamic_time_step!(model::FESeaIceModel, Δt) = dynamic_time_step!(model, model.advection, Δt)

function dynamic_time_step!(model::FESeaIceModel, advection, Δt)
    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ℵ  = model.ice_concentration
    hs = model.snow_thickness
    tracers = model.tracers

    Gⁿ = model.timestepper.Gⁿ

    launch!(arch, grid, :xy, _dynamic_step_tracers!, h, ℵ, h, ℵ, hs, hs, tracers, Gⁿ, Δt)

    return nothing
end

# Thickness and concentration are updated
# We compute hⁿ⁺¹ and ℵⁿ⁺¹ in the same kernel to account for ridging:
# if ℵ > 1, we reset the concentration to 1 and adjust the thickness
# to conserve the total ice volume in the cell.
@kernel function _dynamic_step_tracers!(h, ℵ, hⁿ, ℵⁿ, hs, hsⁿ, tracers, Gⁿ, Δt)
    i, j = @index(Global, NTuple)
    k = 1

    Ghⁿ = Gⁿ.h
    Gℵⁿ = Gⁿ.ℵ

    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁺ = hⁿ[i, j, k] + Δt * Ghⁿ[i, j, k]
        ℵ⁺ = ℵⁿ[i, j, k] + Δt * Gℵⁿ[i, j, k]

        ℵ⁺ = max(zero(ℵ⁺), ℵ⁺) # Concentration cannot be negative, clip it up
        h⁺ = max(zero(h⁺), h⁺) # Thickness cannot be negative, clip it up

        ℵ⁺ = ifelse(h⁺ == 0, zero(ℵ⁺), ℵ⁺) # reset the concentration if there is no sea-ice
        h⁺ = ifelse(ℵ⁺ == 0, zero(h⁺), h⁺) # reset the thickness if there is no sea-ice

        # Ridging and rafting caused by the advection step
        V⁺ = h⁺ * ℵ⁺

        ℵ[i, j, k] = ifelse(ℵ⁺ > 1, one(ℵ⁺), ℵ⁺)
        h[i, j, k] = ifelse(ℵ⁺ > 1, V⁺, h⁺)
    end

    dynamic_step_snow!(i, j, k, hs, hsⁿ, ℵ, Gⁿ, Δt)
end

@inline dynamic_step_snow!(i, j, k, ::Nothing, args...) = nothing

@inline function dynamic_step_snow!(i, j, k, hs, hsⁿ, ℵ, Gⁿ, Δt)
    @inbounds begin
        hs⁺ = hsⁿ[i, j, k] + Δt * Gⁿ.hs[i, j, k]
        hs⁺ = max(zero(hs⁺), hs⁺)
        hs⁺ = ifelse(ℵ[i, j, k] ≤ 0, zero(hs⁺), hs⁺)
        hs[i, j, k] = hs⁺
    end
    return nothing
end
