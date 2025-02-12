using Adapt
using Oceananigans.Fields: ZeroField

# Fallback
@inline implicit_τx_coefficient(i, j, k, grid, stress, clock, fields) = zero(grid)
@inline implicit_τy_coefficient(i, j, k, grid, stress, clock, fields) = zero(grid)

# Fallback
@inline explicit_τx(i, j, k, grid, stress, clock, fields) = zero(grid)
@inline explicit_τy(i, j, k, grid, stress, clock, fields) = zero(grid)

@inline explicit_τx(i, j, k, grid, stress::Number, clock, fields) = stress
@inline explicit_τy(i, j, k, grid, stress::Number, clock, fields) = stress

@inline explicit_τx(i, j, k, grid, stress::AbstractArray, clock, fields) =  @inbounds stress[i, j, k] 
@inline explicit_τy(i, j, k, grid, stress::AbstractArray, clock, fields) =  @inbounds stress[i, j, k] 

# NamedTuple stess (assuming it is `u` and `v`)
@inline implicit_τx_coefficient(i, j, k, grid, stress::NamedTuple, clock, fields) = implicit_τx_coefficient(i, j, k, grid, stress.u, clock, fields)
@inline implicit_τy_coefficient(i, j, k, grid, stress::NamedTuple, clock, fields) = implicit_τy_coefficient(i, j, k, grid, stress.v, clock, fields)

@inline explicit_τx(i, j, k, grid, stress::NamedTuple, clock, fields) = explicit_τx(i, j, k, grid, stress.u, clock, fields)
@inline explicit_τy(i, j, k, grid, stress::NamedTuple, clock, fields) = explicit_τx(i, j, k, grid, stress.v, clock, fields)

#####
##### SemiImplicitStress
#####

struct SemiImplicitStress{U, V, FT}
    uₑ :: U
    vₑ :: V
    ρₑ :: FT
    Cᴰ :: FT
end

"""
    SemiImplicitStress(FT = Float64; 
                       uₑ = ZeroField(FT), 
                       vₑ = ZeroField(FT), 
                       ρₑ = 1026.0, 
                       Cᴰ = 5.5e-3)

A structure representing the semi-implicit stress between the sea ice and an external fluid (either the ocean or the atmosphere),
calculated as
```math
τᵤ = ρₑ Cᴰ sqrt((uₑ - uᵢⁿ)² + (vₑ - vᵢⁿ)²) (uₑ - uᵢⁿ⁺¹)
```
```math
τᵥ = ρₑ Cᴰ sqrt((uₑ - uᵢⁿ)² + (vₑ - vᵢⁿ)²) (vₑ - vᵢⁿ⁺¹)
```

where `uₑ` and `vₑ` are the external velocities, `uᵢⁿ` and `vᵢⁿ` are the sea ice velocities at the current time step,
and `uᵢⁿ⁺¹` and `vᵢⁿ⁺¹` are the sea ice velocities at the next time step.

Arguments
==========
- `FT`: The field type of the velocities (optional, default: Float64).

Keyword Arguments
==================
- `uₑ`: The external x-velocity field.
- `vₑ`: The external y-velocity field.
- `ρₑ`: The density of the external fluid.
- `Cᴰ`: The drag coefficient.
"""
function SemiImplicitStress(FT = Float64; 
                            uₑ = ZeroField(FT), 
                            vₑ = ZeroField(FT), 
                            ρₑ = 1026.0, 
                            Cᴰ = 5.5e-3) 

    return SemiImplicitStress(uₑ, vₑ, convert(FT, ρₑ), convert(FT, Cᴰ))
end

Adapt.adapt_structure(to, τ::SemiImplicitStress) = 
               SemiImplicitStress(Adapt.adapt(to, τ.uₑ), 
                                  Adapt.adapt(to, τ.vₑ), 
                                  τ.ρₑ,
                                  τ.Cᴰ)

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
