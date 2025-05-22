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

# Just passing ocean velocities without mitigation
@inline free_drift_u(i, j, k, grid, f::NamedTuple, clock, model_fields)  = @inbounds f.u[i, j, k] 
@inline free_drift_v(i, j, k, grid, f::NamedTuple, clock, model_fields)  = @inbounds f.v[i, j, k] 

# Passing no velocities
@inline free_drift_u(i, j, k, grid, ::Nothing, clock, model_fields) = zero(grid)
@inline free_drift_v(i, j, k, grid, ::Nothing, clock, model_fields) = zero(grid)

# What if we want to use _only_ the free drift velocities? (not advised)
function time_step_momentum!(model, dynamics::FreeDriftModel, args...)
    
    model_fields = fields(model)
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy, _free_drift_velocities_step!, model, dynamics, model_fields)

    fill_halo_regions!(model.velocities)
    mask_immersed_field_xy!(model.velocities.u, k=size(grid, 3))
    mask_immersed_field_xy!(model.velocities.v, k=size(grid, 3))

    return nothing
end

@kernel function _step_free_drift!(u, v, grid, dynamics, model_fields)
    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    uᶠ = free_drift_u(i, j, kᴺ, grid, dynamics, model_fields.clock, model_fields)
    vᶠ = free_drift_v(i, j, kᴺ, grid, dynamics, model_fields.clock, model_fields)

    @inbounds begin
        u[i, j, 1] = uᶠ
        v[i, j, 1] = vᶠ
    end

    return nothing
end