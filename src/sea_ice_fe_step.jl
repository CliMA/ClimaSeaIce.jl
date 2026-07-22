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
    hc = model.ice_consolidation_thickness

    launch!(arch, grid, :xy, _dynamic_step_tracers!, h, ℵ, h, ℵ, hs, hs, hc, tracers, Gⁿ, Δt)

    return nothing
end

# Thickness and concentration are advanced by flux-form advection of three conserved integrals — the
# intensive content `𝓋 = ℵ·h`, the concentration `ℵ`, and the thickness `h`:
#
#     Gⁿ.𝓋  ≡ ∂𝓋/∂t = -∇·(U·𝓋)     with `𝓋 = ℵ·h`   (the conserved ice content)
#     Gⁿ.ℵ  ≡ ∂ℵ/∂t = -∇·(U·ℵ)
#     Gⁿ.h  ≡ ∂h/∂t = -∇·(U·h)      (advected thickness, used only as a recovery bound)
#
@kernel function _dynamic_step_tracers!(h, ℵ, hⁿ, ℵⁿ, hs, hsⁿ, hc, tracers, Gⁿ, Δt)
    i, j = @index(Global, NTuple)
    k = 1

    G𝓋ⁿ = Gⁿ.𝓋
    Ghⁿ = Gⁿ.h
    Gℵⁿ = Gⁿ.ℵ

    @inbounds begin
        # Advect content `𝓋 = ℵ·h` and concentration ℵ; clip content undershoot and ℵ into [0, 1].
        𝓋ⁿ = hⁿ[i, j, k] * ℵⁿ[i, j, k]
        𝓋⁺ = max(zero(𝓋ⁿ), 𝓋ⁿ + Δt * G𝓋ⁿ[i, j, k])
        ℵ⁺ = min(max(ℵⁿ[i, j, k] + Δt * Gℵⁿ[i, j, k], zero(𝓋ⁿ)), one(𝓋ⁿ))
        hᵗ = hⁿ[i, j, k] + Δt * Ghⁿ[i, j, k]
        
        empty = 𝓋⁺ ≤ zero(𝓋ⁿ)

        # Cap at the advected thickness `hᵗ` (bounds `𝓋/ℵ` as ℵ→0) and floor at the consolidation
        # thickness `hc` (consolidates thin ice, keeps `h > 0`).
        h⁺ = ifelse(empty, zero(𝓋⁺), max(min(𝓋⁺ / ℵ⁺, hᵗ), hc[i, j, k]))
        ℵ⁺ = ifelse(empty, zero(𝓋⁺), 𝓋⁺ / h⁺)

        # Ridging: fold excess concentration (ℵ⁺ > 1 under convergence) into thickness, conserving 𝓋⁺.
        h⁺ = ifelse(ℵ⁺ > 1, h⁺ * ℵ⁺, h⁺)
        ℵ⁺ = min(ℵ⁺, one(ℵ⁺))
        h⁺ = ifelse(iszero(ℵ⁺), zero(h⁺), h⁺)

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
