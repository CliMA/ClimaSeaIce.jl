using Oceananigans.Advection: EnergyConserving, EnstrophyConserving
using Oceananigans.Coriolis: ActiveWeightedEnergyConserving, ActiveWeightedEnstrophyConserving,
                             AbstractRotation, TriadScheme, fᶠᶠᵃ, x_f_cross_U, y_f_cross_U
using Oceananigans.Grids: immersed_peripheral_node, peripheral_node
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Operators: Δx_qᶜᶠᶜ, Δy_qᶠᶜᶜ, Δx⁻¹ᶠᶜᶜ, Δy⁻¹ᶜᶠᶜ, ℑxᶜᵃᵃ, ℑyᵃᶜᵃ

#####
##### Two-dimensional Coriolis acceleration
#####

@inline masked_Δx_qᶜᶠᶜ(i, j, k, grid, q) = Δx_qᶜᶠᶜ(i, j, k, grid, q)
@inline masked_Δy_qᶠᶜᶜ(i, j, k, grid, q) = Δy_qᶠᶜᶜ(i, j, k, grid, q)

@inline function masked_Δx_qᶜᶠᶜ(i, j, k, ibg::ImmersedBoundaryGrid, q)
    active = !immersed_peripheral_node(i, j, k, ibg, Center(), Face(), Center())
    return ifelse(active, Δx_qᶜᶠᶜ(i, j, k, ibg, q), zero(ibg))
end

@inline function masked_Δy_qᶠᶜᶜ(i, j, k, ibg::ImmersedBoundaryGrid, q)
    active = !immersed_peripheral_node(i, j, k, ibg, Face(), Center(), Center())
    return ifelse(active, Δy_qᶠᶜᶜ(i, j, k, ibg, q), zero(ibg))
end

@inline not_peripheral_nodeᶜᶠᶜ(i, j, k, grid) = !peripheral_node(i, j, k, grid, Center(), Face(), Center())
@inline not_peripheral_nodeᶠᶜᶜ(i, j, k, grid) = !peripheral_node(i, j, k, grid, Face(), Center(), Center())

# Coriolis formulations that do not carry one of the schemes below (and that therefore need a
# vertical velocity) keep Oceananigans' own implementation.
@inline x_f_cross_U_2D(i, j, k, grid, coriolis, U) = x_f_cross_U(i, j, k, grid, coriolis, U)
@inline y_f_cross_U_2D(i, j, k, grid, coriolis, U) = y_f_cross_U(i, j, k, grid, coriolis, U)

@inline x_f_cross_U_2D(i, j, k, grid, ::Nothing, U) = zero(grid)
@inline y_f_cross_U_2D(i, j, k, grid, ::Nothing, U) = zero(grid)

#####
##### Enstrophy-conserving scheme
#####

const ESC = AbstractRotation{<:EnstrophyConserving}

@inline x_f_cross_U_2D(i, j, k, grid, coriolis::ESC, U) =
    - ℑyᵃᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis) * ℑxyᶠᶜᵃ(i, j, k, grid, masked_Δx_qᶜᶠᶜ, U.v) * Δx⁻¹ᶠᶜᶜ(i, j, k, grid)

@inline y_f_cross_U_2D(i, j, k, grid, coriolis::ESC, U) =
    + ℑxᶜᵃᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis) * ℑxyᶜᶠᵃ(i, j, k, grid, masked_Δy_qᶠᶜᶜ, U.u) * Δy⁻¹ᶜᶠᶜ(i, j, k, grid)

#####
##### Energy-conserving scheme
#####

const ENC = AbstractRotation{<:EnergyConserving}

@inline f_ℑx_Δx_vᶠᶠᶜ(i, j, k, grid, coriolis, v) = fᶠᶠᵃ(i, j, k, grid, coriolis) * ℑxᶠᵃᵃ(i, j, k, grid, masked_Δx_qᶜᶠᶜ, v)
@inline f_ℑy_Δy_uᶠᶠᶜ(i, j, k, grid, coriolis, u) = fᶠᶠᵃ(i, j, k, grid, coriolis) * ℑyᵃᶠᵃ(i, j, k, grid, masked_Δy_qᶠᶜᶜ, u)

@inline x_f_cross_U_2D(i, j, k, grid, coriolis::ENC, U) =
    - ℑyᵃᶜᵃ(i, j, k, grid, f_ℑx_Δx_vᶠᶠᶜ, coriolis, U.v) * Δx⁻¹ᶠᶜᶜ(i, j, k, grid)

@inline y_f_cross_U_2D(i, j, k, grid, coriolis::ENC, U) =
    + ℑxᶜᵃᵃ(i, j, k, grid, f_ℑy_Δy_uᶠᶠᶜ, coriolis, U.u) * Δy⁻¹ᶜᶠᶜ(i, j, k, grid)

#####
##### Active-weighted schemes (Jamart and Ozer, 1986)
#####

const AESC = AbstractRotation{<:ActiveWeightedEnstrophyConserving}

@inline function x_f_cross_U_2D(i, j, k, grid, coriolis::AESC, U)
    active_nodes = ℑxyᶠᶜᵃ(i, j, k, grid, not_peripheral_nodeᶜᶠᶜ)
    fv = - ℑyᵃᶜᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis) * ℑxyᶠᶜᵃ(i, j, k, grid, masked_Δx_qᶜᶠᶜ, U.v)
    return ifelse(active_nodes == 0, zero(grid), fv / active_nodes) * Δx⁻¹ᶠᶜᶜ(i, j, k, grid)
end

@inline function y_f_cross_U_2D(i, j, k, grid, coriolis::AESC, U)
    active_nodes = ℑxyᶜᶠᵃ(i, j, k, grid, not_peripheral_nodeᶠᶜᶜ)
    fu = ℑxᶜᵃᵃ(i, j, k, grid, fᶠᶠᵃ, coriolis) * ℑxyᶜᶠᵃ(i, j, k, grid, masked_Δy_qᶠᶜᶜ, U.u)
    return ifelse(active_nodes == 0, zero(grid), fu / active_nodes) * Δy⁻¹ᶜᶠᶜ(i, j, k, grid)
end

const AENC = AbstractRotation{<:ActiveWeightedEnergyConserving}

@inline function x_f_cross_U_2D(i, j, k, grid, coriolis::AENC, U)
    active_nodes = ℑxyᶠᶜᵃ(i, j, k, grid, not_peripheral_nodeᶜᶠᶜ)
    fv = - ℑyᵃᶜᵃ(i, j, k, grid, f_ℑx_Δx_vᶠᶠᶜ, coriolis, U.v) * Δx⁻¹ᶠᶜᶜ(i, j, k, grid)
    return ifelse(active_nodes == 0, zero(grid), fv / active_nodes)
end

@inline function y_f_cross_U_2D(i, j, k, grid, coriolis::AENC, U)
    active_nodes = ℑxyᶜᶠᵃ(i, j, k, grid, not_peripheral_nodeᶠᶜᶜ)
    fu = ℑxᶜᵃᵃ(i, j, k, grid, f_ℑy_Δy_uᶠᶠᶜ, coriolis, U.u) * Δy⁻¹ᶜᶠᶜ(i, j, k, grid)
    return ifelse(active_nodes == 0, zero(grid), fu / active_nodes)
end

#####
##### Triad scheme (Arakawa and Lamb, 1981)
#####

@inline 𝒯⁺⁺(i, j, k, grid, coriolis) = fᶠᶠᵃ(i,   j+1, k, grid, coriolis) + fᶠᶠᵃ(i+1, j+1, k, grid, coriolis) + fᶠᶠᵃ(i+1, j,   k, grid, coriolis)
@inline 𝒯⁻⁺(i, j, k, grid, coriolis) = fᶠᶠᵃ(i,   j,   k, grid, coriolis) + fᶠᶠᵃ(i,   j+1, k, grid, coriolis) + fᶠᶠᵃ(i+1, j+1, k, grid, coriolis)
@inline 𝒯⁺⁻(i, j, k, grid, coriolis) = fᶠᶠᵃ(i+1, j+1, k, grid, coriolis) + fᶠᶠᵃ(i+1, j,   k, grid, coriolis) + fᶠᶠᵃ(i,   j,   k, grid, coriolis)
@inline 𝒯⁻⁻(i, j, k, grid, coriolis) = fᶠᶠᵃ(i+1, j,   k, grid, coriolis) + fᶠᶠᵃ(i,   j,   k, grid, coriolis) + fᶠᶠᵃ(i,   j+1, k, grid, coriolis)

const TS = AbstractRotation{<:TriadScheme}

@inline x_f_cross_U_2D(i, j, k, grid, coriolis::TS, U) =
    - Δx⁻¹ᶠᶜᶜ(i, j, k, grid) / 12 * (𝒯⁺⁺(i-1, j, k, grid, coriolis) * masked_Δx_qᶜᶠᶜ(i-1, j+1, k, grid, U.v) +
                                     𝒯⁻⁺(i,   j, k, grid, coriolis) * masked_Δx_qᶜᶠᶜ(i,   j,   k, grid, U.v) +
                                     𝒯⁺⁻(i-1, j, k, grid, coriolis) * masked_Δx_qᶜᶠᶜ(i-1, j,   k, grid, U.v) +
                                     𝒯⁻⁻(i,   j, k, grid, coriolis) * masked_Δx_qᶜᶠᶜ(i,   j+1, k, grid, U.v))

@inline y_f_cross_U_2D(i, j, k, grid, coriolis::TS, U) =
    + Δy⁻¹ᶜᶠᶜ(i, j, k, grid) / 12 * (𝒯⁻⁻(i, j,   k, grid, coriolis) * masked_Δy_qᶠᶜᶜ(i,   j,   k, grid, U.u) +
                                     𝒯⁺⁺(i, j-1, k, grid, coriolis) * masked_Δy_qᶠᶜᶜ(i+1, j-1, k, grid, U.u) +
                                     𝒯⁻⁺(i, j-1, k, grid, coriolis) * masked_Δy_qᶠᶜᶜ(i,   j-1, k, grid, U.u) +
                                     𝒯⁺⁻(i, j,   k, grid, coriolis) * masked_Δy_qᶠᶜᶜ(i+1, j,   k, grid, U.u))
