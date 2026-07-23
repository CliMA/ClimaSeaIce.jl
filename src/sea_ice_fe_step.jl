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

function dynamic_time_step!(model::FESeaIceModel, Δt)
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

# Concentration `ℵ` and content `𝓋 = ℵ·h` are advanced by flux-form advection; the thickness is then
# recovered as `h = 𝓋/ℵ`. The content flux is thickness-weighted (see `advective_thickness_flux_x`) so
# `𝓋 → 0` together with `ℵ`, keeping the recovered thickness bounded.
#
#     Gⁿ.𝓋  ≡ ∂𝓋/∂t = -∇·(U·ℵ·h)   (the conserved ice content)
#     Gⁿ.ℵ  ≡ ∂ℵ/∂t = -∇·(U·ℵ)
#
@kernel function _dynamic_step_tracers!(h, ℵ, hⁿ, ℵⁿ, hs, hsⁿ, tracers, Gⁿ, Δt)
    i, j = @index(Global, NTuple)
    k = 1

    G𝓋ⁿ = Gⁿ.𝓋
    Gℵⁿ = Gⁿ.ℵ

    @inbounds begin
        𝓋ⁿ = hⁿ[i, j, k] * ℵⁿ[i, j, k]
        𝓋⁺ = max(zero(𝓋ⁿ), 𝓋ⁿ + Δt * G𝓋ⁿ[i, j, k])
        ℵ⁺ = max(zero(𝓋ⁿ), ℵⁿ[i, j, k] + Δt * Gℵⁿ[i, j, k])

        empty = 𝓋⁺ ≤ zero(𝓋ⁿ)
        h⁺ = ifelse(ℵ⁺ > 0, 𝓋⁺ / ℵ⁺, zero(𝓋⁺))

        # Ridging: cap concentration at 1 and fold the excess into thickness, conserving 𝓋⁺.
        h⁺ = ifelse(ℵ⁺ > 1, 𝓋⁺, h⁺)
        ℵ⁺ = min(ℵ⁺, one(ℵ⁺))

        h[i, j, k] = ifelse(empty, zero(h⁺), h⁺)
        ℵ[i, j, k] = ifelse(empty, zero(ℵ⁺), ℵ⁺)
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
