using Adapt
using Oceananigans.Fields: ZeroField

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

materialize_stress(stress, grid) = stress

#####
##### Compute stress coefficients
#####

@inline compute_implicit_stress_coefficients!(i, j, k, grid, stress, args...) = nothing

#####
##### SemiImplicitStress
#####

struct SemiImplicitStress{TU, TV, U, V, FT}
    τᵢᵤ :: TU
    τᵢᵥ :: TV
    uₑ  :: U
    vₑ  :: V
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
    return SemiImplicitStress(nothing, nothing, uₑ, vₑ, convert(FT, ρₑ), convert(FT, Cᴰ))
end

function materialize_stress(τ::SemiImplicitStress, grid)
    τᵢᵤ = Field{Face, Center, Nothing}(grid)
    τᵢᵥ = Field{Center, Face, Nothing}(grid)
    return SemiImplicitStress(τᵢᵤ, τᵢᵥ, τ.uₑ, τ.vₑ, τ.ρₑ, τ.Cᴰ)
end

Adapt.adapt_structure(to, τ::SemiImplicitStress) =
               SemiImplicitStress(Adapt.adapt(to, τ.τᵢᵤ),
                                  Adapt.adapt(to, τ.τᵢᵥ),
                                  Adapt.adapt(to, τ.uₑ),
                                  Adapt.adapt(to, τ.vₑ),
                                  τ.ρₑ,
                                  τ.Cᴰ)

function Base.show(io::IO, τ::SemiImplicitStress)
    print(io, "SemiImplicitStress", '\n')
    print(io, "├── τᵢᵤ: ", summary(τ.τᵢᵤ), '\n')
    print(io, "├── τᵢᵥ: ", summary(τ.τᵢᵥ), '\n')
    print(io, "├── uₑ:  ", summary(τ.uₑ), '\n')
    print(io, "├── vₑ:  ", summary(τ.vₑ), '\n')
    print(io, "├── ρₑ:  ", τ.ρₑ, '\n')
    print(io, "└── Cᴰ:  ", τ.Cᴰ)
end

@inline function explicit_τx(i, j, k, grid, τ::SemiImplicitStress, clock, fields)
    uₑ = @inbounds τ.uₑ[i, j, k]
    Δu = @inbounds τ.uₑ[i, j, k] - fields.u[i, j, k]
    Δv = ℑxyᶠᶜᵃ(i, j, k, grid, τ.vₑ) - ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    return τ.ρₑ * τ.Cᴰ * sqrt(Δu^2 + Δv^2) * uₑ
end

@inline function explicit_τy(i, j, k, grid, τ::SemiImplicitStress, clock, fields)
    vₑ = @inbounds τ.vₑ[i, j, k]
    Δv = @inbounds τ.vₑ[i, j, k] - fields.v[i, j, k]
    Δu = ℑxyᶜᶠᵃ(i, j, k, grid, τ.uₑ) - ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)
    return τ.ρₑ * τ.Cᴰ * sqrt(Δu^2 + Δv^2) * vₑ
end

@inline implicit_τx_coefficient(i, j, k, grid, stress::SemiImplicitStress, args...) = @inbounds stress.τᵢᵤ[i, j, k]
@inline implicit_τy_coefficient(i, j, k, grid, stress::SemiImplicitStress, args...) = @inbounds stress.τᵢᵥ[i, j, k]

@inline function compute_implicit_stress_coefficients!(i, j, k, grid, τ::SemiImplicitStress, clock, fields) 
    Δuᶠᶜᶜ = @inbounds τ.uₑ[i, j, k] - fields.u[i, j, k]
    Δvᶠᶜᶜ = ℑxyᶠᶜᵃ(i, j, k, grid, τ.vₑ) - ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)

    Δuᶜᶠᶜ = ℑxyᶜᶠᵃ(i, j, k, grid, τ.uₑ) - ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)
    Δvᶜᶠᶜ = @inbounds τ.vₑ[i, j, k] - fields.v[i, j, k]

    @inbounds τ.τᵢᵤ[i, j, k] = τ.ρₑ * τ.Cᴰ * sqrt(Δuᶠᶜᶜ^2 + Δvᶠᶜᶜ^2)
    @inbounds τ.τᵢᵥ[i, j, k] = τ.ρₑ * τ.Cᴰ * sqrt(Δuᶜᶠᶜ^2 + Δvᶜᶠᶜ^2)

    return nothing
end
