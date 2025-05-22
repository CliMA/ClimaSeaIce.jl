abstract type AbstractFreeDriftDynamics end

"""
    StressBalanceFreeDrift{T, B, C}

A free drift parameterization that computes the free drift velocities as a balance between
top and bottom stresses ``τa ≈ τo`` where we split the stresses into a linear implicit part
and an explicit part  ``τaˣ = τaˣₑ + u * τaˣᵢ`` and ``τoˣ = τoˣₑ + u * τoˣᵢ`` (and similarly for the y-component) such that 
```
uᶠ = (τoˣₑ - τaˣₑ) / (τoˣᵢ - τaˣᵢ)
vᶠ = (τoʸₑ - τaʸₑ) / (τoʸᵢ - τaʸᵢ)
```

Can be used to limit the sea ice velocity when the mass or the concentration are below a certain threshold, or
as a `dynamics` model itself that substitutes the sea ice momentum equation calculation everywhere.
"""
struct StressBalanceFreeDrift{T, B} <: AbstractFreeDriftDynamics
    top_momentum_stess :: T
    bottom_momentum_stress :: B
end

Adapt.adapt_structure(to, s::StressBalanceFreeDrift) = 
    StressBalanceFreeDrift(Adapt.adapt(to, s.top_momentum_stess),
                           Adapt.adapt(to, s.bottom_momentum_stress))
                           
fields(::StressBalanceFreeDrift) = NamedTuple()

@inline function free_drift_u(i, j, k, grid, f::StressBalanceFreeDrift, clock, fields)
    τib = implicit_τx_coefficient(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τit = implicit_τx_coefficient(i, j, k, grid, f.top_momentum_stress, clock, fields)

    τeb = explicit_τx(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τet = explicit_τx(i, j, k, grid, f.top_momentum_stress, clock, fields)

    # Explicit stress
    τe = τeb - τet

    # Implicit stress component 
    τi = τib - τit

    return ifelse(τe == 0, zero(grid), τi / τe)
end

@inline function free_drift_v(i, j, k, grid, f::AbstractFreeDriftDynamics, clock, fields) 
    τib = implicit_τy_coefficient(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τit = implicit_τy_coefficient(i, j, k, grid, f.top_momentum_stress, clock, fields)

    τeb = explicit_τy(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τet = explicit_τy(i, j, k, grid, f.top_momentum_stress, clock, fields)

    # Explicit stress
    τe = τeb - τet

    # Implicit stress component 
    τi = τib - τit

    return ifelse(τe == 0, zero(grid), τi / τe)
end

# Just passing velocities without mitigation
@inline free_drift_u(i, j, k, grid, f::NamedTuple, clock, model_fields)  = @inbounds f.u[i, j, k] 
@inline free_drift_v(i, j, k, grid, f::NamedTuple, clock, model_fields)  = @inbounds f.v[i, j, k] 

# Passing no velocities
@inline free_drift_u(i, j, k, grid, ::Nothing, clock, model_fields) = zero(grid)
@inline free_drift_v(i, j, k, grid, ::Nothing, clock, model_fields) = zero(grid)

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