using Adapt
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: ZeroField, interior
using Oceananigans.Grids: halo_size

# Default no-op implicit stress coefficients for stress types without
# an implicit component.
@inline implicit_τx_coefficient(i, j, k, grid, stress, clock, fields) = zero(grid)
@inline implicit_τy_coefficient(i, j, k, grid, stress, clock, fields) = zero(grid)

# Default no-op explicit stresses for stress types without
# an explicit component.
@inline explicit_τx(i, j, k, grid, stress, clock, fields) = zero(grid)
@inline explicit_τy(i, j, k, grid, stress, clock, fields) = zero(grid)

@inline explicit_τx(i, j, k, grid, stress::Number, clock, fields) = stress
@inline explicit_τy(i, j, k, grid, stress::Number, clock, fields) = stress

@inline explicit_τx(i, j, k, grid, stress::AbstractArray, clock, fields) =  @inbounds stress[i, j, k]
@inline explicit_τy(i, j, k, grid, stress::AbstractArray, clock, fields) =  @inbounds stress[i, j, k]

# NamedTuple stress (assuming it is `u` and `v`)
@inline implicit_τx_coefficient(i, j, k, grid, stress::NamedTuple, clock, fields) = implicit_τx_coefficient(i, j, k, grid, stress.u, clock, fields)
@inline implicit_τy_coefficient(i, j, k, grid, stress::NamedTuple, clock, fields) = implicit_τy_coefficient(i, j, k, grid, stress.v, clock, fields)

@inline explicit_τx(i, j, k, grid, stress::NamedTuple, clock, fields) = explicit_τx(i, j, k, grid, stress.u, clock, fields)
@inline explicit_τy(i, j, k, grid, stress::NamedTuple, clock, fields) = explicit_τy(i, j, k, grid, stress.v, clock, fields)

#####
##### Utility for computing the total stress
#####

@inline x_momentum_stress(i, j, k, grid, stress, clock, fields) =
    @inbounds explicit_τx(i, j, k, grid, stress, clock, fields) - implicit_τx_coefficient(i, j, k, grid, stress, clock, fields) * fields.u[i, j, k]

@inline y_momentum_stress(i, j, k, grid, stress, clock, fields) =
    @inbounds explicit_τy(i, j, k, grid, stress, clock, fields) - implicit_τy_coefficient(i, j, k, grid, stress, clock, fields) * fields.v[i, j, k]

#####
##### Stress materialization
#####

# Whether `field` already lives on `grid`, so it can be read across `grid`'s halo as-is. 
@inline grids_match(field, grid) = field.grid === grid || (field.grid == grid && halo_size(field.grid) == halo_size(grid))

# By default a stress is left untouched (e.g. `nothing`, a `Number`, a `ZeroField`).
extended_external_variable(src, grid) = src

function extended_external_variable(src::Field, grid)
    grids_match(src, grid) && return src
    field = Field{Oceananigans.location(src)...}(grid; boundary_conditions = src.boundary_conditions)
    interior(field) .= interior(src)
    return field
end

# Refresh an extended external velocity from its source, then fill its halo.
function refresh_and_fill_external_velocity!(dst::Field, src)
    interior(dst) .= interior(src)
    fill_halo_regions!(dst)
    return nothing
end

materialize_stress(stress, grid) = extended_external_variable(stress, grid)

function materialize_stress(stress::NamedTuple, grid)
    u = extended_external_variable(stress.u, grid)
    v = extended_external_variable(stress.v, grid)
    return (; u, v)
end

# Fill the external stresses' halos once per time step, before substepping (coupler owns interiors).
update_external_stress!(stress, grid) = nothing

function update_external_stress!(stress::NamedTuple, grid)
    stress.u isa Field && fill_halo_regions!(stress.u)
    stress.v isa Field && fill_halo_regions!(stress.v)
    return nothing
end

#####
##### SemiImplicitStress
#####

struct SemiImplicitStress{U, V, US, VS, FT}
    uₑ  :: U   # external x-velocity read by the kernel (an extended copy after materialization)
    vₑ  :: V   # external y-velocity read by the kernel
    uₑ₀ :: US  # source x-velocity (e.g. live ocean-surface field); copied into uₑ each time step
    vₑ₀ :: VS  # source y-velocity
    ρₑ  :: FT
    Cᴰ  :: FT
end

"""
    SemiImplicitStress(FT = Oceananigans.defaults.FloatType;
                       uₑ = ZeroField(FT),
                       vₑ = ZeroField(FT),
                       ρₑ = 1026.0,
                       Cᴰ = 5.5e-3)

A structure representing the semi-implicit stress between the sea ice and an external fluid
(either the ocean or the atmosphere), calculated as:

```math
\\begin{align*}
τᵤ & = ρₑ Cᴰ \\sqrt{(uₑ - uᵢⁿ)² + (vₑ - vᵢⁿ)²} (uₑ - uᵢⁿ⁺¹) \\\\
τᵥ & = ρₑ Cᴰ \\sqrt{(uₑ - uᵢⁿ)² + (vₑ - vᵢⁿ)²} (vₑ - vᵢⁿ⁺¹)
\\end{align*}
```

where ``uₑ`` and ``vₑ`` are the external velocities, ``uᵢⁿ`` and ``vᵢⁿ`` are the sea ice velocities
at the current time step, and ``uᵢⁿ⁺¹`` and ``vᵢⁿ⁺¹`` are the sea ice velocities at the next time step.

Arguments
==========
- `FT`: The field type of the velocities (optional, default: Oceananigans.defaults.FloatType).

Keyword Arguments
==================
- `uₑ`: The external x-velocity field.
- `vₑ`: The external y-velocity field.
- `ρₑ`: The density of the external fluid.
- `Cᴰ`: The drag coefficient.
"""
function SemiImplicitStress(FT = Oceananigans.defaults.FloatType;
                            uₑ = ZeroField(FT),
                            vₑ = ZeroField(FT),
                            ρₑ = 1026.0,
                            Cᴰ = 5.5e-3)
    return SemiImplicitStress(uₑ, vₑ, uₑ, vₑ, convert(FT, ρₑ), convert(FT, Cᴰ))
end

function materialize_stress(τ::SemiImplicitStress, grid)
    # Extended copies of the external velocities, refreshed from the original source each time step.
    uₑ = extended_external_variable(τ.uₑ₀, grid)
    vₑ = extended_external_variable(τ.vₑ₀, grid)
    return SemiImplicitStress(uₑ, vₑ, τ.uₑ₀, τ.vₑ₀, τ.ρₑ, τ.Cᴰ)
end

# drop source velocities on the device.
Adapt.adapt_structure(to, τ::SemiImplicitStress) =
    SemiImplicitStress(Adapt.adapt(to, τ.uₑ),
                       Adapt.adapt(to, τ.vₑ),
                       nothing,
                       nothing,
                       τ.ρₑ,
                       τ.Cᴰ)

function update_external_stress!(τ::SemiImplicitStress, grid)
    τ.uₑ === τ.uₑ₀ || refresh_and_fill_external_velocity!(τ.uₑ, τ.uₑ₀)
    τ.vₑ === τ.vₑ₀ || refresh_and_fill_external_velocity!(τ.vₑ, τ.vₑ₀)
    return nothing
end

function Base.show(io::IO, τ::SemiImplicitStress)
    print(io, "SemiImplicitStress", '\n')
    print(io, "├── uₑ:  ", summary(τ.uₑ), '\n')
    print(io, "├── vₑ:  ", summary(τ.vₑ), '\n')
    print(io, "├── ρₑ:  ", τ.ρₑ, '\n')
    print(io, "└── Cᴰ:  ", τ.Cᴰ)
end

@inline function x_momentum_stress(i, j, k, grid, τ::SemiImplicitStress, clock, fields)
    uₑ = @inbounds τ.uₑ[i, j, k]
    Δu = @inbounds uₑ - fields.u[i, j, k]
    Δv = ℑxyᶠᶜᵃ(i, j, k, grid, τ.vₑ) - ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    return τ.ρₑ * τ.Cᴰ * sqrt(Δu^2 + Δv^2) * Δu
end

@inline function y_momentum_stress(i, j, k, grid, τ::SemiImplicitStress, clock, fields)
    vₑ = @inbounds τ.vₑ[i, j, k]
    Δv = @inbounds vₑ - fields.v[i, j, k]
    Δu = ℑxyᶜᶠᵃ(i, j, k, grid, τ.uₑ) - ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)
    return τ.ρₑ * τ.Cᴰ * sqrt(Δu^2 + Δv^2) * Δv
end

@inline function explicit_τx(i, j, k, grid, τ::SemiImplicitStress, clock, fields)
    uₑ = @inbounds τ.uₑ[i, j, k]
    Δu = @inbounds uₑ - fields.u[i, j, k]
    Δv = ℑxyᶠᶜᵃ(i, j, k, grid, τ.vₑ) - ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    return τ.ρₑ * τ.Cᴰ * sqrt(Δu^2 + Δv^2) * uₑ
end

@inline function explicit_τy(i, j, k, grid, τ::SemiImplicitStress, clock, fields)
    vₑ = @inbounds τ.vₑ[i, j, k]
    Δv = @inbounds vₑ - fields.v[i, j, k]
    Δu = ℑxyᶜᶠᵃ(i, j, k, grid, τ.uₑ) - ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)
    return τ.ρₑ * τ.Cᴰ * sqrt(Δu^2 + Δv^2) * vₑ
end

# Computed on the fly from the current velocities so that, inside the alternating substep,
# v's drag sees the just-updated u (and vice versa).
@inline function implicit_τx_coefficient(i, j, k, grid, τ::SemiImplicitStress, clock, fields)
    Δu = @inbounds τ.uₑ[i, j, k] - fields.u[i, j, k]
    Δv = ℑxyᶠᶜᵃ(i, j, k, grid, τ.vₑ) - ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    return τ.ρₑ * τ.Cᴰ * sqrt(Δu^2 + Δv^2)
end

@inline function implicit_τy_coefficient(i, j, k, grid, τ::SemiImplicitStress, clock, fields)
    Δu = ℑxyᶜᶠᵃ(i, j, k, grid, τ.uₑ) - ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)
    Δv = @inbounds τ.vₑ[i, j, k] - fields.v[i, j, k]
    return τ.ρₑ * τ.Cᴰ * sqrt(Δu^2 + Δv^2)
end
