using Adapt

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

"""
    struct SemiImplicitOceanSeaIceStress{U, V, FT}

A structure representing the semi-implicit stress between the ocean and sea ice, 
calculated as
```math
τᵤ = ρₒ Cᴰ [(uₒ - uᵢⁿ)² + (vₒ - vᵢⁿ)²] (uₒ - uᵢⁿ⁺¹)
τᵥ = ρₒ Cᴰ [(uₒ - uᵢⁿ)² + (vₒ - vᵢⁿ)²] (vₒ - vᵢⁿ⁺¹)
```

Fields
======
- `u`: The zonal (x-direction) component of the ocean velocity.
- `v`: The meridional (y-direction) component of the ocean velocity.
- `ρₒ`: The density of the ocean.
- `Cᴰ`: The drag coefficient between the ocean and sea ice.
"""
struct SemiImplicitOceanSeaIceStress{U, V, FT}
    u  :: U
    v  :: V
    ρₒ :: FT
    Cᴰ :: FT
end

Adapt.adapt_structure(to, τ::SemiImplicitOceanSeaIceStress) = 
    SemiImplicitOceanSeaIceStress(Adapt.adapt(to, τ.u), 
                                  Adapt.adapt(to, τ.v), 
                                  τ.ρₒ,
                                  τ.Cᴰ)

@inline function explicit_τx(i, j, k, grid, τ::SemiImplicitOceanSeaIceStress, clock, fields) 
    uₒ = @inbounds τ.u[i, j, k]
    Δu = @inbounds fields.u[i, j, k] - τ.u[i, j, k]
    Δv = ℑxyᶠᶜᵃ(i, j, k, grid, τ.v) - ℑxyᶠᶜᵃ(i, j, k, grid, fields.v) 
    return τ.ρₒCᴰ * sqrt(Δu^2 + Δv^2) * uₒ
end

@inline function explicit_τy(i, j, k, grid, τ::SemiImplicitOceanSeaIceStress, clock, fields) 
    vₒ = @inbounds τ.v[i, j, k]
    Δu = ℑxyᶜᶠᵃ(i, j, k, grid, τ.u) - ℑxyᶜᶠᵃ(i, j, k, grid, fields.u) 
    Δv = @inbounds fields.v[i, j, k] - τ.v[i, j, k] 
    return τ.ρₒCᴰ * sqrt(Δu^2 + Δv^2) * vₒ
end

@inline function implicit_τx_coefficient(i, j, k, grid, τ::SemiImplicitOceanSeaIceStress, clock, fields) 
    Δu = @inbounds fields.u[i, j, k] - τ.u[i, j, k]
    Δv = ℑxyᶠᶜᵃ(i, j, k, grid, τ.v) - ℑxyᶠᶜᵃ(i, j, k, grid, fields.v) 
    return τ.ρₒCᴰ * sqrt(Δu^2 + Δv^2)
end

@inline function implicit_τy_coefficient(i, j, k, grid, τ::SemiImplicitOceanSeaIceStress, clock, fields) 
    Δu = ℑxyᶜᶠᵃ(i, j, k, grid, τ.u) - ℑxyᶜᶠᵃ(i, j, k, grid, fields.u) 
    Δv = @inbounds fields.v[i, j, k] - τ.v[i, j, k] 
    return τ.ρₒCᴰ * sqrt(Δu^2 + Δv^2)
end
