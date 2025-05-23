abstract type AbstractFreeDriftDynamics end

"""
    StressBalanceFreeDrift{T, B}

A free drift parameterization that computes the free drift velocities as a balance between top and bottom stresses ``τa ≈ τo``.

In case the only one of the stresses is a `SemiImplicitStress`, the model will compute the free drift velocity exactly 
assuming that the other stress does not depend on the sea ice velocity. Otherwise, the model will compute the free drift 
velocity approximately using the previous time step velocities to compute the nonlinear terms.

Can be used to limit the sea ice velocity when the mass or the concentration are below a certain threshold, or
as a `dynamics` model itself that substitutes the sea ice momentum equation calculation everywhere.
"""
struct StressBalanceFreeDrift{T, B} <: AbstractFreeDriftDynamics
    top_momentum_stress :: T
    bottom_momentum_stress :: B
end

Adapt.adapt_structure(to, s::StressBalanceFreeDrift) = 
    StressBalanceFreeDrift(Adapt.adapt(to, s.top_momentum_stess),
                           Adapt.adapt(to, s.bottom_momentum_stress))

fields(::StressBalanceFreeDrift) = NamedTuple()

# Stress balance when either the top or the bottom stresses do not depend on ice velocity
# In this case we have a simplified form of the free drift velocity
# Otherwise, to avoid a nonlinear solve, we assume the stress is only lineary dependent on the velocity at time-step
# n+1 and use the ice velocities at time-step n to compute the nonlinear term. 
# Note that this is the same formulation we use to solve for stresses in the `SeaIceMomentumEquation` dynamics.
const TISB = StressBalanceFreeDrift{<:Union{AbstractArray, NamedTuple}, <:SemiImplicitStress}
const BISB = StressBalanceFreeDrift{<:SemiImplicitStress, <:Union{AbstractArray, NamedTuple}}

# Stress balance when only the bottom stress is ice-velocity dependent:
# Then: 𝒰ᵢ = 𝒰ᴮ - τᵀ / sqrt(Cᴮ * ||τᵀ||)
@inline function free_drift_u(i, j, k, grid, f::TISB, clock, fields) 
    τxᵀ = x_momentum_stress(i, j, k, grid, f.top_momentum_stress, clock, fields)
    τyᵀ = y_momentum_stress(i, j, k, grid, f.top_momentum_stress, clock, fields)
    τᵀ  = sqrt(τxᵀ^2 + τyᵀ^2)

    τᴮ = f.bottom_momentum_stress
    uᴮ = @inbounds τᴮ.uₑ[i, j, k]
    Cᴮ = τᴮ.ρₑ * τᴮ.Cᴰ

    return uᴮ - ifelse(τᵀ == 0, τᵀ, τxᵀ / sqrt(Cᴮ * τᵀ))
end

@inline function free_drift_v(i, j, k, grid, f::TISB, clock, fields) 
    τxᵀ = x_momentum_stress(i, j, k, grid, f.top_momentum_stress, clock, fields)
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
    τyᴮ = y_momentum_stress(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τᴮ  = sqrt(τxᴮ^2 + τyᴮ^2)

    τᵀ = f.top_momentum_stess
    uᵀ = @inbounds τᵀ.uₑ[i, j, k]
    Cᵀ = τᵀ.ρₑ * τᵀ.Cᴰ

    return uᵀ - ifelse(τᴮ == 0, τᴮ, τxᴮ / sqrt(Cᵀ * τᴮ))
end

@inline function free_drift_v(i, j, k, grid, f::BISB, clock, fields) 
    τxᴮ = x_momentum_stress(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τyᴮ = y_momentum_stress(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τᴮ  = sqrt(τxᴮ^2 + τyᴮ^2)

    τᵀ = f.top_momentum_stess
    vᵀ = @inbounds τᵀ.vₑ[i, j, k]
    Cᵀ = τᵀ.ρₑ * τᵀ.Cᴰ

    return vᵀ - ifelse(τᴮ == 0, τᴮ, τyᴮ / sqrt(Cᵀ * τᴮ))
end

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
