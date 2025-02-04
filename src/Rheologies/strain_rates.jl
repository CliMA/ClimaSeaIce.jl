using Oceananigans.ImmersedBoundaries: IBG

@inline strain_rate_xx(i, j, k, grid, u, v) = δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, u) / Azᶜᶜᶜ(i, j, k, grid)
@inline strain_rate_yy(i, j, k, grid, u, v) = δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, v) / Azᶜᶜᶜ(i, j, k, grid)
@inline strain_rate_xy(i, j, k, grid, u, v) = (δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶠᶜ, v) + δyᵃᶠᵃ(i, j, k, grid, Δx_qᶠᶜᶜ, u)) / Azᶠᶠᶜ(i, j, k, grid) / 2

@inline strain_rate_xx(i, j, k, grid::IBG, u, v) = ∂xᴮᶜᶜᶜ(i, j, k, grid, u)
@inline strain_rate_yy(i, j, k, grid::IBG, u, v) = ∂yᴮᶜᶜᶜ(i, j, k, grid, v)
@inline strain_rate_xy(i, j, k, grid::IBG, u, v) = (∂xᴮᶠᶠᶜ(i, j, k, grid, u) + ∂yᴮᶠᶠᶜ(i, j, k, grid, v)) / 2

@inline function strain_rate_xy_centered(i, j, k, grid, scheme, u, v)
   δv = δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, ℑxyᶠᶜᵃ, v) 
   δu = δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, ℑxyᶜᶠᵃ, u)
   return (δv + δu) / Azᶜᶜᶜ(i, j, k, grid) / 2
end

# Hardcode No-slip boundary conditions on immersed boundaries?
# TODO: Fix this mess!! 
# Find a way not to hard-code. We need to pass the immersed_boundary_conditions of
# the velocities to the kernels
@inline function ∂xᴮᶜᶜᶜ(i, j, k, grid, u)
    i1 = inactive_node(i,   j, k, grid, f, c, c)
    i2 = inactive_node(i+1, j, k, grid, f, c, c) 
    Az = Azᶜᶜᶜ(i, j, k, grid)

    u1 = @inbounds u[i,   j, k] * Δyᶠᶜᶜ(i,   j, k, grid)
    u2 = @inbounds u[i+1, j, k] * Δyᶠᶜᶜ(i+1, j, k, grid)

    return ifelse(i1, 2u2 / Az, ifelse(i2, - 2u1 / Az, (u2 - u1) / Az))
end

@inline function ∂yᴮᶜᶜᶜ(i, j, k, grid, v)
    j1 = inactive_node(i, j,   k, grid, c, f, c)
    j2 = inactive_node(i, j+1, k, grid, c, f, c) 
    Az = Azᶜᶜᶜ(i, j, k, grid)

    v1 = @inbounds v[i, j,   k] * Δxᶜᶠᶜ(i,   j, k, grid)
    v2 = @inbounds v[i, j+1, k] * Δxᶜᶠᶜ(i+1, j, k, grid)

    return ifelse(j1, 2v2 / Az, ifelse(j2, 2v1 / Az, (v2 - v1) / Az))
end

@inline function ∂xᴮᶠᶠᶜ(i, j, k, grid, v)
    i1 = inactive_node(i-1, j, k, grid, c, f, c)
    i2 = inactive_node(i,   j, k, grid, c, f, c) 
    Az = Azᶠᶠᶜ(i, j, k, grid)

    v1 = @inbounds v[i-1, j, k] * Δyᶜᶠᶜ(i,   j, k, grid)
    v2 = @inbounds v[i,   j, k] * Δyᶜᶠᶜ(i+1, j, k, grid)

    return ifelse(i1, 2v2 / Az, ifelse(i2, - 2v1 / Az, (v2 - v1) / Az))
end

@inline function ∂yᴮᶠᶠᶜ(i, j, k, grid, u)
    j1 = inactive_node(i, j-1, k, grid, f, c, c)
    j2 = inactive_node(i, j,   k, grid, f, c, c) 
    Az = Azᶠᶠᶜ(i, j, k, grid)

    u1 = @inbounds u[i, j-1, k] * Δxᶠᶜᶜ(i,   j, k, grid)
    u2 = @inbounds u[i, j,   k] * Δxᶠᶜᶜ(i+1, j, k, grid)

    return ifelse(j1, 2u2 / Az, ifelse(j2, 2u1 / Az, (u2 - u1) / Az))
end
