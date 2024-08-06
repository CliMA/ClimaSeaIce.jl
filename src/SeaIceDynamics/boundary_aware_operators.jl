using Oceananigans.Grids: AbstractUnderlyingGrid
using Oceananigans.ImmersedBoundaries
using Oceananigans.ImmersedBoundaries: inactive_node
using Oceananigans.Operators

# These operators are crucial for sea-ice rheologies and explicit solvers that operate 
# on a C-grid because they require a looooot of interpolations. For this reason it is
# better to have a the dynamics solved on an E-grid (next PR)
# TODO: Move them to Oceananigans?

const AUG = AbstractUnderlyingGrid
const IBG = ImmersedBoundaryGrid

const c = Center()
const f = Face()

@inline ı(i, j, k, grid, F::Function, args...) = F(i, j, k, grid, args...)
@inline ı(i, j, k, grid, ϕ)                    = @inbounds ϕ[i, j, k]

# Defining Interpolation operators for the immersed boundaries
@inline conditional_ℑx_f(LY, LZ, i, j, k, grid, ℑx, args...) = ifelse(inactive_node(i, j, k, grid, c, LY, LZ), ı(i-1, j, k, grid, args...), ifelse(inactive_node(i-1, j, k, grid, c, LY, LZ), ı(i, j, k, grid, args...), ℑx(i, j, k, grid, args...)))
@inline conditional_ℑx_c(LY, LZ, i, j, k, grid, ℑx, args...) = ifelse(inactive_node(i, j, k, grid, f, LY, LZ), ı(i+1, j, k, grid, args...), ifelse(inactive_node(i+1, j, k, grid, f, LY, LZ), ı(i, j, k, grid, args...), ℑx(i, j, k, grid, args...)))
@inline conditional_ℑy_f(LX, LZ, i, j, k, grid, ℑy, args...) = ifelse(inactive_node(i, j, k, grid, LX, c, LZ), ı(i, j-1, k, grid, args...), ifelse(inactive_node(i, j-1, k, grid, LX, c, LZ), ı(i, j, k, grid, args...), ℑy(i, j, k, grid, args...)))
@inline conditional_ℑy_c(LX, LZ, i, j, k, grid, ℑy, args...) = ifelse(inactive_node(i, j, k, grid, LX, f, LZ), ı(i, j+1, k, grid, args...), ifelse(inactive_node(i, j+1, k, grid, LX, f, LZ), ı(i, j, k, grid, args...), ℑy(i, j, k, grid, args...)))

@inline ℑxᴮᶜᶜᶜ(i, j, k, grid, args...) = conditional_ℑx_c(c, c, i, j, k, grid, ℑxᶜᵃᵃ, args...)
@inline ℑxᴮᶠᶜᶜ(i, j, k, grid, args...) = conditional_ℑx_f(c, c, i, j, k, grid, ℑxᶠᵃᵃ, args...)
@inline ℑyᴮᶜᶜᶜ(i, j, k, grid, args...) = conditional_ℑy_c(c, c, i, j, k, grid, ℑyᵃᶜᵃ, args...)
@inline ℑyᴮᶜᶠᶜ(i, j, k, grid, args...) = conditional_ℑy_f(c, c, i, j, k, grid, ℑyᵃᶠᵃ, args...)

@inline ℑxᴮᶜᶠᶜ(i, j, k, grid, args...) = conditional_ℑx_c(f, c, i, j, k, grid, ℑxᶜᵃᵃ, args...)
@inline ℑxᴮᶠᶠᶜ(i, j, k, grid, args...) = conditional_ℑx_f(f, c, i, j, k, grid, ℑxᶠᵃᵃ, args...)
@inline ℑyᴮᶠᶜᶜ(i, j, k, grid, args...) = conditional_ℑy_c(f, c, i, j, k, grid, ℑyᵃᶜᵃ, args...)
@inline ℑyᴮᶠᶠᶜ(i, j, k, grid, args...) = conditional_ℑy_f(f, c, i, j, k, grid, ℑyᵃᶠᵃ, args...)

@inline ℑxyᴮᶜᶜᶜ(i, j, k, grid, args...) = ℑxᴮᶜᶜᶜ(i, j, k, grid, ℑyᴮᶜᶜᶜ, args...)
@inline ℑxyᴮᶠᶜᶜ(i, j, k, grid, args...) = ℑxᴮᶠᶜᶜ(i, j, k, grid, ℑyᴮᶠᶜᶜ, args...)
@inline ℑxyᴮᶜᶠᶜ(i, j, k, grid, args...) = ℑxᴮᶜᶠᶜ(i, j, k, grid, ℑyᴮᶜᶠᶜ, args...)
@inline ℑxyᴮᶠᶠᶜ(i, j, k, grid, args...) = ℑxᴮᶠᶠᶜ(i, j, k, grid, ℑyᴮᶠᶠᶜ, args...)

using Oceananigans.Grids: AbstractGrid

const AG{X, Y, Z} = AbstractGrid{<:Any, X, Y, Z} where {X, Y, Z}
const P = Periodic
const B = Bounded
const R = RightConnected
const L = LeftConnected

# The functions `η★` `U★` and `V★` represent the value of free surface, barotropic zonal and meridional velocity at time step m+1/2
@inline δxᶠᵃᵃ_c(i, j, k, grid::AG, η) = δxᶠᵃᵃ(i, j, k, grid, η)
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG, η) = δyᵃᶠᵃ(i, j, k, grid, η)
@inline δxᶜᵃᵃ_U(i, j, k, grid::AG, U) = δxᶜᵃᵃ(i, j, k, grid, U)
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG, V) = δyᵃᶜᵃ(i, j, k, grid, V)

@inline δxᶠᵃᵃ_c(i, j, k, grid::AG{<:P},        η) = @inbounds ifelse(i == 1, η[1, j, k] - η[grid.Nx, j, k], δxᶠᵃᵃ(i, j, k, grid, η))
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG{<:Any, <:P}, η) = @inbounds ifelse(j == 1, η[i, 1, k] - η[i, grid.Ny, k], δyᵃᶠᵃ(i, j, k, grid, η))

@inline δxᶜᵃᵃ_U(i, j, k, grid::AG{<:P},        U) = ifelse(i == grid.Nx, U[1, j, k] - U[grid.Nx, j, k], δxᶜᵃᵃ(i, j, k, grid, U))
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG{<:Any, <:P}, V) = ifelse(j == grid.Ny, V[i, 1, k] - V[i, grid.Ny, k], δyᵃᶜᵃ(i, j, k, grid, V))

# Enforce NoFlux conditions for `η★`
@inline δxᶠᵃᵃ_c(i, j, k, grid::AG{<:B},        η) = ifelse(i == 1, zero(grid), δxᶠᵃᵃ(i, j, k, grid, η))
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG{<:Any, <:B}, η) = ifelse(j == 1, zero(grid), δyᵃᶠᵃ(i, j, k, grid, η))
@inline δxᶠᵃᵃ_c(i, j, k, grid::AG{<:R},        η) = ifelse(i == 1, zero(grid), δxᶠᵃᵃ(i, j, k, grid, η))
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG{<:Any, <:R}, η) = ifelse(j == 1, zero(grid), δyᵃᶠᵃ(i, j, k, grid, η))

# Enforce Impenetrability conditions for `U★` and `V★`
@inline δxᶜᵃᵃ_U(i, j, k, grid::AG{<:B},        U) = ifelse(i == grid.Nx, - U[i, j, k], ifelse(i == 1, U[2, j, k], δxᶜᵃᵃ(i, j, k, grid, U)))
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG{<:Any, <:B}, V) = ifelse(j == grid.Ny, - V[i, j, k], ifelse(j == 1, V[i, 2, k], δyᵃᶜᵃ(i, j, k, grid, V)))

@inline δxᶜᵃᵃ_U(i, j, k, grid::AG{<:L},        U) = ifelse(i == grid.Nx, - U[i, j, k], δxᶜᵃᵃ(i, j, k, grid, U))
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG{<:Any, <:L}, V) = ifelse(j == grid.Ny, - V[i, j, k], δyᵃᶜᵃ(i, j, k, grid, V))

@inline δxᶜᵃᵃ_U(i, j, k, grid::AG{<:R},        U) = ifelse(i == 1, U[2, j, k], δxᶜᵃᵃ(i, j, k, grid, U))
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG{<:Any, <:R}, V) = ifelse(j == 1, V[i, 2, k], δyᵃᶜᵃ(i, j, k, grid, V))

# Derivative Operators
@inline ∂xᶠᵃᵃ_c(i, j, k, grid, η) = δxᶠᵃᵃ_c(i, j, k, grid, η) / Δxᶠᶜᶜ(i, j, k, grid)
@inline ∂yᵃᶠᵃ_c(i, j, k, grid, η) = δyᵃᶠᵃ_c(i, j, k, grid, η) / Δyᶜᶠᶜ(i, j, k, grid)

# Derivative Operators
@inline ∂xᶜᵃᵃ_U(i, j, k, grid, U) = δxᶜᵃᵃ_U(i, j, k, grid, U) / Δxᶜᶜᶜ(i, j, k, grid)
@inline ∂yᵃᶜᵃ_V(i, j, k, grid, V) = δyᵃᶜᵃ_V(i, j, k, grid, V) / Δyᶜᶜᶜ(i, j, k, grid)

using Oceananigans.ImmersedBoundaries: conditional_δx_f, 
                                       conditional_δx_c,
                                       conditional_δy_f,
                                       conditional_δy_c,
                                       translate_loc

const c = Center()
const f = Face()

@inline ∂xᶠᶠᶜ_c(i, j, k, grid, args...) = ∂xᶠᵃᵃ_c(i, j, k, grid, args...)
@inline ∂yᶠᶠᶜ_c(i, j, k, grid, args...) = ∂yᵃᶠᵃ_c(i, j, k, grid, args...)
@inline ∂xᶜᶠᶜ_U(i, j, k, grid, args...) = ∂xᶜᵃᵃ_U(i, j, k, grid, args...)
@inline ∂yᶠᶜᶜ_V(i, j, k, grid, args...) = ∂yᵃᶜᵃ_V(i, j, k, grid, args...)

@inline ∂xᶠᶠᶜ_c(i, j, k, grid::IBG, args...) = conditional_δx_f(f, c, i, j, k, grid, ∂xᶠᵃᵃ_c, args...)
@inline ∂yᶠᶠᶜ_c(i, j, k, grid::IBG, args...) = conditional_δy_f(f, c, i, j, k, grid, ∂yᵃᶠᵃ_c, args...)
@inline ∂xᶜᶠᶜ_U(i, j, k, grid::IBG, args...) = conditional_δx_c(f, c, i, j, k, grid, ∂xᶜᵃᵃ_U, args...)
@inline ∂yᶠᶜᶜ_V(i, j, k, grid::IBG, args...) = conditional_δy_c(f, c, i, j, k, grid, ∂yᵃᶜᵃ_V, args...)

@inline ∂xᶠᶜᶜ_c(i, j, k, grid, args...) = ∂xᶠᵃᵃ_c(i, j, k, grid, args...)
@inline ∂yᶜᶠᶜ_c(i, j, k, grid, args...) = ∂yᵃᶠᵃ_c(i, j, k, grid, args...)
@inline ∂xᶜᶜᶜ_U(i, j, k, grid, args...) = ∂xᶜᵃᵃ_U(i, j, k, grid, args...)
@inline ∂yᶜᶜᶜ_V(i, j, k, grid, args...) = ∂yᵃᶜᵃ_V(i, j, k, grid, args...)

@inline ∂xᶠᶜᶜ_c(i, j, k, grid::IBG, args...) = conditional_δx_f(c, c, i, j, k, grid, ∂xᶠᵃᵃ_c, args...)
@inline ∂yᶜᶠᶜ_c(i, j, k, grid::IBG, args...) = conditional_δy_f(c, c, i, j, k, grid, ∂yᵃᶠᵃ_c, args...)
@inline ∂xᶜᶜᶜ_U(i, j, k, grid::IBG, args...) = conditional_δx_c(c, c, i, j, k, grid, ∂xᶜᵃᵃ_U, args...)
@inline ∂yᶜᶜᶜ_V(i, j, k, grid::IBG, args...) = conditional_δy_c(c, c, i, j, k, grid, ∂yᵃᶜᵃ_V, args...)