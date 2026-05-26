using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, 
                                       ImmersedBoundaryCondition

using Oceananigans.BoundaryConditions: FBC, getbc
using Oceananigans.Advection:conditional_flux_ccc, conditional_flux_ffc
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

@inline _ice_stress_ux(i, j, k, ibg::IBG, args...) = conditional_flux_ccc(i, j, k, ibg, zero(ibg), ice_stress_ux(i, j, k, ibg, args...))
@inline _ice_stress_uy(i, j, k, ibg::IBG, args...) = conditional_flux_ffc(i, j, k, ibg, zero(ibg), ice_stress_uy(i, j, k, ibg, args...))
@inline _ice_stress_vx(i, j, k, ibg::IBG, args...) = conditional_flux_ffc(i, j, k, ibg, zero(ibg), ice_stress_vx(i, j, k, ibg, args...))
@inline _ice_stress_vy(i, j, k, ibg::IBG, args...) = conditional_flux_ccc(i, j, k, ibg, zero(ibg), ice_stress_vy(i, j, k, ibg, args...))

#####
##### Stress divergence
#####

@inline Δy²_qᶜᶜᶜ(i, j, k, grid, q, args...) = Δyᶜᶜᶜ(i, j, k, grid)^2 * q(i, j, k, grid, args...)
@inline Δx²_qᶠᶠᶜ(i, j, k, grid, q, args...) = Δxᶠᶠᶜ(i, j, k, grid)^2 * q(i, j, k, grid, args...)
@inline Δy²_qᶠᶠᶜ(i, j, k, grid, q, args...) = Δyᶠᶠᶜ(i, j, k, grid)^2 * q(i, j, k, grid, args...)
@inline Δx²_qᶜᶜᶜ(i, j, k, grid, q, args...) = Δxᶜᶜᶜ(i, j, k, grid)^2 * q(i, j, k, grid, args...)

# Stress invariants at their native locations
@inline σI( i, j, k, grid, args...) = _ice_stress_ux(i, j, k, grid, args...) + _ice_stress_vy(i, j, k, grid, args...)
@inline σII(i, j, k, grid, args...) = _ice_stress_ux(i, j, k, grid, args...) - _ice_stress_vy(i, j, k, grid, args...)

@inline function ∂ⱼ_σ₁ⱼ(i, j, k, grid, rheology, clock, fields)
    δ = Δyᶠᶜᶜ(i, j, k, grid) * δxᶠᵃᵃ(i, j, k, grid, σI, rheology, clock, fields) / 2
    T = δxᶠᵃᵃ(i, j, k, grid, Δy²_qᶜᶜᶜ, σII, rheology, clock, fields) / Δyᶠᶜᶜ(i, j, k, grid) / 2
    S = δyᵃᶜᵃ(i, j, k, grid, Δx²_qᶠᶠᶜ, _ice_stress_uy, rheology, clock, fields) / Δxᶠᶜᶜ(i, j, k, grid)
    return (δ + T + S) / Azᶠᶜᶜ(i, j, k, grid)
end

@inline function ∂ⱼ_σ₂ⱼ(i, j, k, grid, rheology, clock, fields)
    δ  =   Δxᶜᶠᶜ(i, j, k, grid) * δyᵃᶠᵃ(i, j, k, grid, σI, rheology, clock, fields) / 2
    T  = - δyᵃᶠᵃ(i, j, k, grid, Δx²_qᶜᶜᶜ, σII, rheology, clock, fields) / Δxᶜᶠᶜ(i, j, k, grid) / 2
    S =    δxᶜᵃᵃ(i, j, k, grid, Δy²_qᶠᶠᶜ, _ice_stress_vx, rheology, clock, fields) / Δyᶜᶠᶜ(i, j, k, grid)
    return (δ + T + S) / Azᶜᶠᶜ(i, j, k, grid)
end

#####
##### Immersed Stress divergence
#####

@inline immersed_∂ⱼ_σ₁ⱼ(i, j, k, grid, args...) = zero(grid)
@inline immersed_∂ⱼ_σ₂ⱼ(i, j, k, grid, args...) = zero(grid)

@inline flip(::Type{Center}) = Face
@inline flip(::Type{Face})   = Center
@inline flip(::Center) = Face()
@inline flip(::Face)   = Center()

@inline function immersed_∂ⱼ_σ₁ⱼ(i, j, k, ibg::IBG, u_bc::IBC, rheology, clock, fields)
    # Fetch fluxes across immersed boundary
    q̃ᵂ =  ib_ice_stress_ux_west(i, j, k, ibg, u_bc.west,  rheology, clock, fields)
    q̃ᴱ =  ib_ice_stress_ux_east(i, j, k, ibg, u_bc.east,  rheology, clock, fields)
    q̃ˢ = ib_ice_stress_uy_south(i, j, k, ibg, u_bc.south, rheology, clock, fields)
    q̃ᴺ = ib_ice_stress_uy_north(i, j, k, ibg, u_bc.north, rheology, clock, fields)

    iᵂ, jˢ, _ = map(index_left,  (i, j, k), (f, c, c)) # Broadcast instead of map causes inference failure
    iᴱ, jᴺ, _ = map(index_right, (i, j, k), (f, c, c))

    # Impose i) immersed fluxes if we're on an immersed boundary or ii) zero otherwise.
    qᵂ = conditional_flux_ccc(iᵂ, j, k, ibg, q̃ᵂ, zero(ibg)) * Axᶜᶜᶜ(iᵂ, j, k, ibg)
    qᴱ = conditional_flux_ccc(iᴱ, j, k, ibg, q̃ᴱ, zero(ibg)) * Axᶜᶜᶜ(iᴱ, j, k, ibg)
    qˢ = conditional_flux_ffc(i, jˢ, k, ibg, q̃ˢ, zero(ibg)) * Ayᶠᶠᶜ(i, jˢ, k, ibg)
    qᴺ = conditional_flux_ffc(i, jᴺ, k, ibg, q̃ᴺ, zero(ibg)) * Ayᶠᶠᶜ(i, jᴺ, k, ibg)

    return (qᴱ - qᵂ + qᴺ - qˢ) / Vᶠᶜᶜ(i, j, k, ibg)
end

@inline function immersed_∂ⱼ_σ₂ⱼ(i, j, k, ibg::IBG, v_bc::IBC, rheology, clock, fields)
    # Fetch fluxes across immersed boundary
    q̃ᵂ =  ib_ice_stress_vx_west(i, j, k, ibg, v_bc.west,  rheology, clock, fields)
    q̃ᴱ =  ib_ice_stress_vx_east(i, j, k, ibg, v_bc.east,  rheology, clock, fields)
    q̃ˢ = ib_ice_stress_vy_south(i, j, k, ibg, v_bc.south, rheology, clock, fields)
    q̃ᴺ = ib_ice_stress_vy_north(i, j, k, ibg, v_bc.north, rheology, clock, fields)

    iᵂ, jˢ, _ = map(index_left,  (i, j, k), (c, f, c)) # Broadcast instead of map causes inference failure
    iᴱ, jᴺ, _ = map(index_right, (i, j, k), (c, f, c))

    # Impose i) immersed fluxes if we're on an immersed boundary or ii) zero otherwise.
    qᵂ = conditional_flux_ffc(iᵂ, j, k, ibg, q̃ᵂ, zero(ibg)) * Axᶠᶠᶜ(iᵂ, j, k, ibg)
    qᴱ = conditional_flux_ffc(iᴱ, j, k, ibg, q̃ᴱ, zero(ibg)) * Axᶠᶠᶜ(iᴱ, j, k, ibg)
    qˢ = conditional_flux_ccc(i, jˢ, k, ibg, q̃ˢ, zero(ibg)) * Ayᶜᶜᶜ(i, jˢ, k, ibg)
    qᴺ = conditional_flux_ccc(i, jᴺ, k, ibg, q̃ᴺ, zero(ibg)) * Ayᶜᶜᶜ(i, jᴺ, k, ibg)

    return (qᴱ - qᵂ + qᴺ - qˢ) / Vᶜᶠᶜ(i, j, k, ibg)
end

# TODO: Implement immersed fluxes (0 for the moment)
@inline  ib_ice_stress_ux_west(i, j, k, ibg, args...) = zero(ibg)
@inline  ib_ice_stress_ux_east(i, j, k, ibg, args...) = zero(ibg)
@inline  ib_ice_stress_vx_west(i, j, k, ibg, args...) = zero(ibg)
@inline  ib_ice_stress_vx_east(i, j, k, ibg, args...) = zero(ibg)
@inline ib_ice_stress_uy_south(i, j, k, ibg, args...) = zero(ibg)
@inline ib_ice_stress_uy_north(i, j, k, ibg, args...) = zero(ibg)
@inline ib_ice_stress_vy_south(i, j, k, ibg, args...) = zero(ibg)
@inline ib_ice_stress_vy_north(i, j, k, ibg, args...) = zero(ibg)

# Only supporting FluxBoundaryConditions for now. The user-supplied `bc` follows the
# Oceananigans convention (flux into the cell along the inward-pointing normal), so we
# negate to recover the stress σ = -F before plugging into the +∇·σ formula.
@inline  ib_ice_stress_ux_west(i, j, k, ibg, bc::FBC, rheology, args...) = - getbc(bc, i, j, k, ibg, args...)
@inline  ib_ice_stress_ux_east(i, j, k, ibg, bc::FBC, rheology, args...) = + getbc(bc, i, j, k, ibg, args...)
@inline  ib_ice_stress_vx_west(i, j, k, ibg, bc::FBC, rheology, args...) = - getbc(bc, i, j, k, ibg, args...)
@inline  ib_ice_stress_vx_east(i, j, k, ibg, bc::FBC, rheology, args...) = + getbc(bc, i, j, k, ibg, args...)
@inline ib_ice_stress_uy_south(i, j, k, ibg, bc::FBC, rheology, args...) = - getbc(bc, i, j, k, ibg, args...)
@inline ib_ice_stress_uy_north(i, j, k, ibg, bc::FBC, rheology, args...) = + getbc(bc, i, j, k, ibg, args...)
@inline ib_ice_stress_vy_south(i, j, k, ibg, bc::FBC, rheology, args...) = - getbc(bc, i, j, k, ibg, args...)
@inline ib_ice_stress_vy_north(i, j, k, ibg, bc::FBC, rheology, args...) = + getbc(bc, i, j, k, ibg, args...)
