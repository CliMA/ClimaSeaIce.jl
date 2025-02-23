using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, 
                                       ImmersedBoundaryCondition,
                                       flip, 
                                       conditional_flux, 
                                       index_left, 
                                       index_right

using Oceananigans.Advection:conditional_flux_ccc, conditional_flux_ffc

const IBG = ImmersedBoundaryGrid
const IBC = ImmersedBoundaryCondition

# Utilities
@inline _ice_stress_ux(args...) = ice_stress_ux(args...)
@inline _ice_stress_uy(args...) = ice_stress_uy(args...)
@inline _ice_stress_vx(args...) = ice_stress_vx(args...)
@inline _ice_stress_vy(args...) = ice_stress_vy(args...)

@inline _ice_stress_ux(i, j, k, ibg::IBG, args...) = conditional_flux_ccc(i, j, k, ibg, zero(grid), ice_stress_ux(i, j, k, ibg, args...))
@inline _ice_stress_uy(i, j, k, ibg::IBG, args...) = conditional_flux_ffc(i, j, k, ibg, zero(grid), ice_stress_uy(i, j, k, ibg, args...))
@inline _ice_stress_vx(i, j, k, ibg::IBG, args...) = conditional_flux_ffc(i, j, k, ibg, zero(grid), ice_stress_vx(i, j, k, ibg, args...))
@inline _ice_stress_vy(i, j, k, ibg::IBG, args...) = conditional_flux_ccc(i, j, k, ibg, zero(grid), ice_stress_vy(i, j, k, ibg, args...))

# Stress divergence
@inline function ∂ⱼ_σ₁ⱼ(i, j, k, grid, rheology, clock, fields)
    return 1 / Azᶠᶜᶜ(i, j, k, grid) * (δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶜᶜ, _ice_stress_ux, rheology, clock, fields) +
                                       δyᵃᶜᵃ(i, j, k, grid, Δx_qᶠᶠᶜ, _ice_stress_uy, rheology, clock, fields))
end

@inline function ∂ⱼ_σ₂ⱼ(i, j, k, grid, rheology, clock, fields)
    return 1 / Azᶠᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶠᶜ, _ice_stress_vx, rheology, clock, fields) +
                                       δyᵃᶠᵃ(i, j, k, grid, Δx_qᶜᶜᶜ, _ice_stress_vy, rheology, clock, fields))
end

# Immersed Stress divergence
@inline immersed_∂ⱼ_σ₁ⱼ(i, j, k, grid, args...) = zero(grid)
@inline immersed_∂ⱼ_σ₂ⱼ(i, j, k, grid, args...) = zero(grid)

@inline immersed_∂ⱼ_σ₁ⱼ(i, j, k, ibg::IBG, ::ZFBC, args...) = zero(ibg)
@inline immersed_∂ⱼ_σ₂ⱼ(i, j, k, ibg::IBG, ::ZFBC, args...) = zero(ibg)

# Value boundary condition implementation
@inline function immersed_∂ⱼ_σ₁ⱼ(i, j, k, ibg::IBG, u_bc::IBC, rheology, clock, fields)
    # Fetch fluxes across immersed boundary
    q̃ᵂ = ib_ice_stress_ux(i,   j,   k, ibg, u_bc.west,  rheology, clock, fields)
    q̃ᴱ = ib_ice_stress_ux(i+1, j,   k, ibg, u_bc.east,  rheology, clock, fields)
    q̃ˢ = ib_ice_stress_uy(i,   j-1, k, ibg, u_bc.south, rheology, clock, fields)
    q̃ᴺ = ib_ice_stress_uy(i,   j,   k, ibg, u_bc.north, rheology, clock, fields)

    iᵂ, jˢ, kᴮ = map(index_left,  (i, j, k), loc) # Broadcast instead of map causes inference failure
    iᴱ, jᴺ, kᵀ = map(index_right, (i, j, k), loc)
    LX, LY, LZ = loc

    # Impose i) immersed fluxes if we're on an immersed boundary or ii) zero otherwise.
    qᵂ = conditional_flux(iᵂ, j, k, ibg, flip(LX), LY, LZ, q̃ᵂ, zero(ibg))
    qᴱ = conditional_flux(iᴱ, j, k, ibg, flip(LX), LY, LZ, q̃ᴱ, zero(ibg))
    qˢ = conditional_flux(i, jˢ, k, ibg, LX, flip(LY), LZ, q̃ˢ, zero(ibg))
    qᴺ = conditional_flux(i, jᴺ, k, ibg, LX, flip(LY), LZ, q̃ᴺ, zero(ibg))

    return div(i, j, k, ibg, loc, qᵂ, qᴱ, qˢ, qᴺ, qᴮ, qᵀ)
end

# Value boundary condition implementation
@inline function immersed_∂ⱼ_σ₂ⱼ(i, j, k, ibg::IBG, v_bc::IBC, rheology, clock, fields)
    # Fetch fluxes across immersed boundary
    q̃ᵂ = ib_ice_stress_vx(i-1, j,   k, ibg, v_bc.west,  rheology, clock, fields)
    q̃ᴱ = ib_ice_stress_vx(i,   j,   k, ibg, v_bc.east,  rheology, clock, fields)
    q̃ˢ = ib_ice_stress_vy(i,   j,   k, ibg, v_bc.south, rheology, clock, fields)
    q̃ᴺ = ib_ice_stress_vy(i,   j+1, k, ibg, v_bc.north, rheology, clock, fields)

    iᵂ, jˢ, kᴮ = map(index_left,  (i, j, k), loc) # Broadcast instead of map causes inference failure
    iᴱ, jᴺ, kᵀ = map(index_right, (i, j, k), loc)
    LX, LY, LZ = loc

    # Impose i) immersed fluxes if we're on an immersed boundary or ii) zero otherwise.
    qᵂ = conditional_flux(iᵂ, j, k, ibg, flip(LX), LY, LZ, q̃ᵂ, zero(ibg))
    qᴱ = conditional_flux(iᴱ, j, k, ibg, flip(LX), LY, LZ, q̃ᴱ, zero(ibg))
    qˢ = conditional_flux(i, jˢ, k, ibg, LX, flip(LY), LZ, q̃ˢ, zero(ibg))
    qᴺ = conditional_flux(i, jᴺ, k, ibg, LX, flip(LY), LZ, q̃ᴺ, zero(ibg))

    return div(i, j, k, ibg, loc, qᵂ, qᴱ, qˢ, qᴺ, qᴮ, qᵀ)
end

# TODO: Implement immersed fluxes
@inline ib_ice_stress_ux(i, j, k, grid, args...) = zero(grid)
@inline ib_ice_stress_vx(i, j, k, grid, args...) = zero(grid)
@inline ib_ice_stress_uy(i, j, k, grid, args...) = zero(grid)
@inline ib_ice_stress_vy(i, j, k, grid, args...) = zero(grid)