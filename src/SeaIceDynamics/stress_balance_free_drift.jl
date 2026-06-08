abstract type AbstractFreeDriftDynamics end

struct StressBalanceFreeDrift{T, B} <: AbstractFreeDriftDynamics
    top_momentum_stress :: T
    bottom_momentum_stress :: B
end

"""
    StressBalanceFreeDrift(; top_momentum_stress = nothing,
                             bottom_momentum_stress = nothing)

A free drift parameterization that computes the free drift velocities as a balance between top and bottom stresses ``τₐ ≈ τₒ``.

The only supported configuration is when either the `top_momentum_stress` or the `bottom_momentum_stress` are a
`SemiImplicitStress`. The model will compute the free drift velocity exactly assuming that the other stress does
not depend on the sea ice velocity.

Can be used to limit the sea ice velocity when the mass or the concentration are below a certain threshold, or
as a `dynamics` model itself that substitutes the sea ice momentum equation calculation everywhere.
"""
function StressBalanceFreeDrift(; top_momentum_stress = nothing,
                                  bottom_momentum_stress = nothing)

    if top_momentum_stress isa SemiImplicitStress
        if bottom_momentum_stress isa SemiImplicitStress
            throw(ArgumentError("`StressBalanceFreeDrift` supports a `SemiImplicitStress` only for the `top_momentum_stress` or the `bottom_momentum_stress`, not both"))
        end
    else
        if !(bottom_momentum_stress isa SemiImplicitStress)
            throw(ArgumentError("`StressBalanceFreeDrift` requires using a `SemiImplicitStress` for either the `top_momentum_stress` or the `bottom_momentum_stress`"))
        end
    end

    return StressBalanceFreeDrift(top_momentum_stress, bottom_momentum_stress)
end

Adapt.adapt_structure(to, s::StressBalanceFreeDrift) =
    StressBalanceFreeDrift(Adapt.adapt(to, s.top_momentum_stress),
                           Adapt.adapt(to, s.bottom_momentum_stress))

# Repoint a free drift at already-materialized external stresses (see `materialize_solver`).
materialize_free_drift(free_drift, top_momentum_stress, bottom_momentum_stress) = free_drift
materialize_free_drift(::StressBalanceFreeDrift, top_momentum_stress, bottom_momentum_stress) = StressBalanceFreeDrift(top_momentum_stress, bottom_momentum_stress)

Oceananigans.fields(::StressBalanceFreeDrift) = NamedTuple()

# Stress balance when exactly one of the top or bottom stresses is
# ice-velocity dependent. In this case we have a simplified closed-form
# expression for the free-drift velocity.
#
# If both stresses are velocity dependent, the free-drift velocity would
# require a nonlinear solve, so we define an explicit mixed-case dispatch
# that throws an informative error.
const TISB = StressBalanceFreeDrift{<:Any, <:SemiImplicitStress}
const BISB = StressBalanceFreeDrift{<:SemiImplicitStress, <:Any}
const TBISB = StressBalanceFreeDrift{<:SemiImplicitStress, <:SemiImplicitStress}

# Stress balance when only the bottom stress is ice-velocity dependent:
# Then: 𝒰ᵢ = 𝒰ᴮ - τᵀ / sqrt(Cᴮ * ||τᵀ||)
@inline function free_drift_u(i, j, k, grid, f::TISB, clock, fields)
    τxᵀ = x_momentum_stress(i, j, k, grid, f.top_momentum_stress, clock, fields)
    τyᵀ = ℑxyᶠᶜᵃ(i, j, k, grid, y_momentum_stress, f.top_momentum_stress, clock, fields)
    τᵀ  = sqrt(τxᵀ^2 + τyᵀ^2)

    τᴮ = f.bottom_momentum_stress
    uᴮ = @inbounds τᴮ.uₑ[i, j, k]
    Cᴮ = τᴮ.ρₑ * τᴮ.Cᴰ

    return uᴮ - ifelse(τᵀ == 0, τᵀ, τxᵀ / sqrt(Cᴮ * τᵀ))
end

@inline function free_drift_v(i, j, k, grid, f::TISB, clock, fields)
    τxᵀ = ℑxyᶜᶠᵃ(i, j, k, grid, x_momentum_stress, f.top_momentum_stress, clock, fields)
    τyᵀ = y_momentum_stress(i, j, k, grid, f.top_momentum_stress, clock, fields)
    τᵀ  = sqrt(τxᵀ^2 + τyᵀ^2)

    τᴮ = f.bottom_momentum_stress
    vᴮ = @inbounds τᴮ.vₑ[i, j, k]
    Cᴮ = τᴮ.ρₑ * τᴮ.Cᴰ

    return vᴮ - ifelse(τᵀ == 0, τᵀ, τyᵀ / sqrt(Cᴮ * τᵀ))
end

# Stress balance when only the top stress is ice-velocity dependent:
# Then: 𝒰ᵢ = 𝒰ᵀ - τᴮ / sqrt(Cᵀ * ||τᴮ||)
@inline function free_drift_u(i, j, k, grid, f::BISB, clock, fields)
    τxᴮ = x_momentum_stress(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τyᴮ = ℑxyᶠᶜᵃ(i, j, k, grid, y_momentum_stress, f.bottom_momentum_stress, clock, fields)
    τᴮ  = sqrt(τxᴮ^2 + τyᴮ^2)

    τᵀ = f.top_momentum_stress
    uᵀ = @inbounds τᵀ.uₑ[i, j, k]
    Cᵀ = τᵀ.ρₑ * τᵀ.Cᴰ

    return uᵀ - ifelse(τᴮ == 0, τᴮ, τxᴮ / sqrt(Cᵀ * τᴮ))
end

@inline function free_drift_v(i, j, k, grid, f::BISB, clock, fields)
    τxᴮ = ℑxyᶜᶠᵃ(i, j, k, grid, x_momentum_stress, f.bottom_momentum_stress, clock, fields)
    τyᴮ = y_momentum_stress(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τᴮ  = sqrt(τxᴮ^2 + τyᴮ^2)

    τᵀ = f.top_momentum_stress
    vᵀ = @inbounds τᵀ.vₑ[i, j, k]
    Cᵀ = τᵀ.ρₑ * τᵀ.Cᴰ

    return vᵀ - ifelse(τᴮ == 0, τᴮ, τyᴮ / sqrt(Cᵀ * τᴮ))
end

@noinline function unsupported_mixed_stress_balance_free_drift()
    throw(ArgumentError("`StressBalanceFreeDrift` does not support both `top_momentum_stress` and `bottom_momentum_stress` being `SemiImplicitStress`. This mixed `TBISB` case requires a nonlinear solve that is not implemented."))
end

@inline free_drift_u(i, j, k, grid, ::TBISB, clock, fields) = unsupported_mixed_stress_balance_free_drift()
@inline free_drift_v(i, j, k, grid, ::TBISB, clock, fields) = unsupported_mixed_stress_balance_free_drift()

const NoFreeDrift = StressBalanceFreeDrift{<:Nothing, <:Nothing}

@inline free_drift_u(i, j, k, grid, ::NoFreeDrift, clock, fields) = zero(grid)
@inline free_drift_v(i, j, k, grid, ::NoFreeDrift, clock, fields) = zero(grid)

# Fallbacks for a given velocity field.
@inline free_drift_u(i, j, k, grid, f::NamedTuple, clock, fields)  = @inbounds f.u[i, j, k]
@inline free_drift_v(i, j, k, grid, f::NamedTuple, clock, fields)  = @inbounds f.v[i, j, k]

# Passing no velocities
@inline free_drift_u(i, j, k, grid, ::Nothing, clock, fields) = zero(grid)
@inline free_drift_v(i, j, k, grid, ::Nothing, clock, fields) = zero(grid)

# What if we want to use _only_ the free drift velocities? (not advised)
function time_step_momentum!(model, dynamics::AbstractFreeDriftDynamics, args...)

    model_fields = fields(model)
    clock = model.clock
    grid  = model.grid
    arch  = architecture(grid)
    u, v  = model.velocities

    launch!(arch, grid, :xy, _free_drift_velocity_step!, u, v, grid, dynamics, clock, model_fields)

    return nothing
end

@kernel function _free_drift_velocity_step!(u, v, grid, dynamics, clock, fields)
    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    @inbounds u[i, j, 1] = free_drift_u(i, j, kᴺ, grid, dynamics, clock, fields)
    @inbounds v[i, j, 1] = free_drift_v(i, j, kᴺ, grid, dynamics, clock, fields)
end
