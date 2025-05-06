using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, 
                                       ImmersedBoundaryCondition

using Oceananigans.Advection:conditional_flux_ccc, conditional_flux_ffc
using Oceananigans.Operators
using Oceananigans.Operators: index_left, index_right

const IBG = ImmersedBoundaryGrid
const IBC = ImmersedBoundaryCondition

const c = Center()
const f = Face()

#####
##### Ice internal stresses
#####

@inline _ice_stress_ux(args...) = ice_stress_ux(args...)
@inline _ice_stress_uy(args...) = ice_stress_uy(args...)
@inline _ice_stress_vx(args...) = ice_stress_vx(args...)
@inline _ice_stress_vy(args...) = ice_stress_vy(args...)

#####
##### Stress divergence
#####

@inline function ∂ⱼ_σ₁ⱼ(i, j, k, grid, rheology, clock, fields)
    return Az⁻¹ᶠᶠᶜ(i, j, k, grid) * (δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶠᶜ, _ice_stress_ux, rheology, clock, fields) +
                                     δyᵃᶠᵃ(i, j, k, grid, Δx_qᶠᶜᶜ, _ice_stress_uy, rheology, clock, fields))
end

@inline function ∂ⱼ_σ₂ⱼ(i, j, k, grid, rheology, clock, fields)
    return Az⁻¹ᶠᶠᶜ(i, j, k, grid) * (δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶠᶜ, _ice_stress_vx, rheology, clock, fields) +
                                     δyᵃᶠᵃ(i, j, k, grid, Δx_qᶠᶜᶜ, _ice_stress_vy, rheology, clock, fields))
end