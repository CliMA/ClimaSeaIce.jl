using Oceananigans.Units: Time
using Oceananigans.Fields: flattened_unique_values, ZeroField
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!
using Oceananigans.ImmersedBoundaries: mask_immersed_field_xy!

using ClimaSeaIce.SeaIceDynamics: time_step_momentum!
using ClimaSeaIce.SeaIceThermodynamics: thermodynamic_time_step!

const FESeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:Any, <:ForwardEulerTimeStepper}

# We separate the thermodynamic step from the advection (dynamic) step.
# The thermodynamic step is column physics and is performed all at once.
function time_step!(model::FESeaIceModel, Δt; kwargs...)
    
    # Be paranoid and update state at iteration 0
    model.clock.iteration == 0 && update_state!(model)

    # Compute advective tendencies and update 
    # advected tracers
    compute_tendencies!(model, Δt)
    dynamic_time_step!(model, Δt)

    # Perform the thermodynamic step
    thermodynamic_time_step!(model, model.ice_thermodynamics, model.snow_thermodynamics, Δt)

    # This is an implicit (or split-explicit) step to advance momentum.
    time_step_momentum!(model, model.dynamics, Δt)

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

    Gⁿ   = model.timestepper.Gⁿ
    ℵmin = model.concentration_floor

    launch!(arch, grid, :xy, _dynamic_step_tracers!, h, ℵ, h, ℵ, hs, hs, ℵmin, tracers, Gⁿ, Δt)

    return nothing
end

# Thickness and concentration are updated using a flux-form advection of the
# intensive ice content per area `𝓋 = ℵ·h`
#
#     Gⁿ.h  ≡ ∂𝓋/∂t = -∇·(U·𝓋)    with 𝓋 = ℵ·h
#     Gⁿ.ℵ  ≡ ∂ℵ/∂t = -∇·(U·ℵ)
#
# `𝓋⁺ = 𝓋ⁿ + Δt · Gⁿ.h` is updated directly, then `h = 𝓋/ℵ` is recovered with
# a small-ℵ guard (`concentration_floor ≈ 1e-10`). Ridging is automatic: if
# `ℵ⁺ > 1` we cap `ℵ = 1` and `h = 𝓋/1 = 𝓋`, conserving the content the
# advection step produced.
@kernel function _dynamic_step_tracers!(h, ℵ, hⁿ, ℵⁿ, hs, hsⁿ, ℵmin, tracers, Gⁿ, Δt)
    i, j = @index(Global, NTuple)
    k = 1

    Ghⁿ = Gⁿ.h
    Gℵⁿ = Gⁿ.ℵ

    @inbounds begin
        # Advect ice content per area `𝓋 = ℵ·h` (units m) and concentration ℵ. 
        𝓋ⁿ = hⁿ[i, j, k] * ℵⁿ[i, j, k]
        𝓋⁺ = 𝓋ⁿ + Δt * Ghⁿ[i, j, k]
        ℵ⁺ = ℵⁿ[i, j, k] + Δt * Gℵⁿ[i, j, k]

        # Clip undershoots.
        𝓋⁺ = max(zero(𝓋⁺), 𝓋⁺)
        ℵ⁺ = max(zero(ℵ⁺), ℵ⁺)

        # Ridging: cap concentration at 1
        ℵ⁺ = min(ℵ⁺, one(ℵ⁺))

        # h recovery `h = 𝓋 / ℵ`, guarded by `model.concentration_floor`.
        # Below the floor both 𝓋 and ℵ are zapped to zero (mass loss limited to O(ℵmin · h)).
        active = ℵ⁺ > ℵmin
        h⁺     = ifelse(active, 𝓋⁺ / ℵ⁺, zero(𝓋⁺))
        ℵ⁺     = ifelse(active, ℵ⁺,      zero(ℵ⁺))

        ℵ[i, j, k] = ℵ⁺
        h[i, j, k] = h⁺
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
