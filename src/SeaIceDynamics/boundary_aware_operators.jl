using Oceananigans.Grids: AbstractUnderlyingGrid
using Oceananigans.ImmersedBoundaries
using Oceananigans.ImmersedBoundaries: inactive_node
using Oceananigans.Operators

struct Slip end
struct NoSlip end

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
@inline ℑx_f(LY, LZ, i, j, k, grid, bc, ℑx, args...) = ifelse(inactive_node(i, j, k, grid, c, LY, LZ), ı(i-1, j, k, grid, args...), ifelse(inactive_node(i-1, j, k, grid, c, LY, LZ), ı(i, j, k, grid, args...), ℑx(i, j, k, grid, args...)))
@inline ℑx_c(LY, LZ, i, j, k, grid, bc, ℑx, args...) = ifelse(inactive_node(i, j, k, grid, f, LY, LZ), ı(i+1, j, k, grid, args...), ifelse(inactive_node(i+1, j, k, grid, f, LY, LZ), ı(i, j, k, grid, args...), ℑx(i, j, k, grid, args...)))
@inline ℑy_f(LX, LZ, i, j, k, grid, bc, ℑy, args...) = ifelse(inactive_node(i, j, k, grid, LX, c, LZ), ı(i, j-1, k, grid, args...), ifelse(inactive_node(i, j-1, k, grid, LX, c, LZ), ı(i, j, k, grid, args...), ℑy(i, j, k, grid, args...)))
@inline ℑy_c(LX, LZ, i, j, k, grid, bc, ℑy, args...) = ifelse(inactive_node(i, j, k, grid, LX, f, LZ), ı(i, j+1, k, grid, args...), ifelse(inactive_node(i, j+1, k, grid, LX, f, LZ), ı(i, j, k, grid, args...), ℑy(i, j, k, grid, args...)))

@inline ℑx_f(LY, LZ, i, j, k, grid, ::NoSlip, ℑx, args...) = ifelse(inactive_node(i, j, k, grid, c, LY, LZ) | inactive_node(i-1, j, k, grid, c, LY, LZ), zero(grid), ℑx(i, j, k, grid, args...))
@inline ℑx_c(LY, LZ, i, j, k, grid, ::NoSlip, ℑx, args...) = ifelse(inactive_node(i, j, k, grid, f, LY, LZ) | inactive_node(i+1, j, k, grid, f, LY, LZ), zero(grid), ℑx(i, j, k, grid, args...))
@inline ℑy_f(LX, LZ, i, j, k, grid, ::NoSlip, ℑy, args...) = ifelse(inactive_node(i, j, k, grid, LX, c, LZ) | inactive_node(i, j-1, k, grid, LX, c, LZ), zero(grid), ℑy(i, j, k, grid, args...))
@inline ℑy_c(LX, LZ, i, j, k, grid, ::NoSlip, ℑy, args...) = ifelse(inactive_node(i, j, k, grid, LX, f, LZ) | inactive_node(i, j+1, k, grid, LX, f, LZ), zero(grid), ℑy(i, j, k, grid, args...))

@inline ℑxᴮᶜᶜᶜ(i, j, k, grid, bc, args...) = ℑx_c(c, c, i, j, k, grid, bc, ℑxᶜᵃᵃ, args...)
@inline ℑxᴮᶠᶜᶜ(i, j, k, grid, bc, args...) = ℑx_f(c, c, i, j, k, grid, bc, ℑxᶠᵃᵃ, args...)
@inline ℑyᴮᶜᶜᶜ(i, j, k, grid, bc, args...) = ℑy_c(c, c, i, j, k, grid, bc, ℑyᵃᶜᵃ, args...)
@inline ℑyᴮᶜᶠᶜ(i, j, k, grid, bc, args...) = ℑy_f(c, c, i, j, k, grid, bc, ℑyᵃᶠᵃ, args...)

@inline ℑxᴮᶜᶠᶜ(i, j, k, grid, bc, args...) = ℑx_c(f, c, i, j, k, grid, bc, ℑxᶜᵃᵃ, args...)
@inline ℑxᴮᶠᶠᶜ(i, j, k, grid, bc, args...) = ℑx_f(f, c, i, j, k, grid, bc, ℑxᶠᵃᵃ, args...)
@inline ℑyᴮᶠᶜᶜ(i, j, k, grid, bc, args...) = ℑy_c(f, c, i, j, k, grid, bc, ℑyᵃᶜᵃ, args...)
@inline ℑyᴮᶠᶠᶜ(i, j, k, grid, bc, args...) = ℑy_f(f, c, i, j, k, grid, bc, ℑyᵃᶠᵃ, args...)

@inline ℑxyᴮᶜᶜᶜ(i, j, k, grid, bc, args...) = ℑxᴮᶜᶜᶜ(i, j, k, grid, bc, ℑyᴮᶠᶜᶜ, bc, args...)
@inline ℑxyᴮᶠᶜᶜ(i, j, k, grid, bc, args...) = ℑxᴮᶠᶜᶜ(i, j, k, grid, bc, ℑyᴮᶜᶜᶜ, bc, args...)
@inline ℑxyᴮᶜᶠᶜ(i, j, k, grid, bc, args...) = ℑxᴮᶜᶠᶜ(i, j, k, grid, bc, ℑyᴮᶠᶠᶜ, bc, args...)
@inline ℑxyᴮᶠᶠᶜ(i, j, k, grid, bc, args...) = ℑxᴮᶠᶠᶜ(i, j, k, grid, bc, ℑyᴮᶜᶠᶜ, bc, args...)

using Oceananigans.Grids: AbstractGrid

const AG{X, Y, Z} = AbstractGrid{<:Any, X, Y, Z} where {X, Y, Z}
const P = Periodic
const B = Bounded
const R = RightConnected
const L = LeftConnected

# The functions `η★` `U★` and `V★` represent the value of free surface, barotropic zonal and meridional velocity at time step m+1/2
@inline δxᶠᵃᵃ_c(i, j, k, grid::AG, bc, η) = δxᶠᵃᵃ(i, j, k, grid, η)
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG, bc, η) = δyᵃᶠᵃ(i, j, k, grid, η)
@inline δxᶜᵃᵃ_U(i, j, k, grid::AG, bc, U) = δxᶜᵃᵃ(i, j, k, grid, U)
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG, bc, V) = δyᵃᶜᵃ(i, j, k, grid, V)

@inline δxᶠᵃᵃ_c(i, j, k, grid::AG{<:P},        bc, η) = @inbounds ifelse(i == 1, η[1, j, k] - η[grid.Nx, j, k], δxᶠᵃᵃ(i, j, k, grid, η))
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG{<:Any, <:P}, bc, η) = @inbounds ifelse(j == 1, η[i, 1, k] - η[i, grid.Ny, k], δyᵃᶠᵃ(i, j, k, grid, η))

@inline δxᶜᵃᵃ_U(i, j, k, grid::AG{<:P},        bc, U) = ifelse(i == grid.Nx, U[1, j, k] - U[grid.Nx, j, k], δxᶜᵃᵃ(i, j, k, grid, U))
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG{<:Any, <:P}, bc, V) = ifelse(j == grid.Ny, V[i, 1, k] - V[i, grid.Ny, k], δyᵃᶜᵃ(i, j, k, grid, V))

# Enforce NoFlux conditions for `η★`
@inline δxᶠᵃᵃ_c(i, j, k, grid::AG{<:B},        bc, η) = ifelse(i == 1, zero(grid), δxᶠᵃᵃ(i, j, k, grid, η))
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG{<:Any, <:B}, bc, η) = ifelse(j == 1, zero(grid), δyᵃᶠᵃ(i, j, k, grid, η))
@inline δxᶠᵃᵃ_c(i, j, k, grid::AG{<:R},        bc, η) = ifelse(i == 1, zero(grid), δxᶠᵃᵃ(i, j, k, grid, η))
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG{<:Any, <:R}, bc, η) = ifelse(j == 1, zero(grid), δyᵃᶠᵃ(i, j, k, grid, η))

# Enforce NoSlip conditions for `η★`
@inline δxᶠᵃᵃ_c(i, j, k, grid::AG{<:B},        ::NoSlip, η) = @inbounds ifelse(i == 1, 2η[i, j, k], δxᶠᵃᵃ(i, j, k, grid, η))
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG{<:Any, <:B}, ::NoSlip, η) = @inbounds ifelse(j == 1, 2η[i, j, k], δyᵃᶠᵃ(i, j, k, grid, η))
@inline δxᶠᵃᵃ_c(i, j, k, grid::AG{<:R},        ::NoSlip, η) = @inbounds ifelse(i == 1, 2η[i, j, k], δxᶠᵃᵃ(i, j, k, grid, η))
@inline δyᵃᶠᵃ_c(i, j, k, grid::AG{<:Any, <:R}, ::NoSlip, η) = @inbounds ifelse(j == 1, 2η[i, j, k], δyᵃᶠᵃ(i, j, k, grid, η))

# Enforce Impenetrability conditions for `U★` and `V★`
@inline δxᶜᵃᵃ_U(i, j, k, grid::AG{<:B},        bc, U) = ifelse(i == grid.Nx, - U[i, j, k], ifelse(i == 1, U[2, j, k], δxᶜᵃᵃ(i, j, k, grid, U)))
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG{<:Any, <:B}, bc, V) = ifelse(j == grid.Ny, - V[i, j, k], ifelse(j == 1, V[i, 2, k], δyᵃᶜᵃ(i, j, k, grid, V)))

@inline δxᶜᵃᵃ_U(i, j, k, grid::AG{<:L},        bc, U) = ifelse(i == grid.Nx, - U[i, j, k], δxᶜᵃᵃ(i, j, k, grid, U))
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG{<:Any, <:L}, bc, V) = ifelse(j == grid.Ny, - V[i, j, k], δyᵃᶜᵃ(i, j, k, grid, V))

@inline δxᶜᵃᵃ_U(i, j, k, grid::AG{<:R},        bc, U) = ifelse(i == 1, U[2, j, k], δxᶜᵃᵃ(i, j, k, grid, U))
@inline δyᵃᶜᵃ_V(i, j, k, grid::AG{<:Any, <:R}, bc, V) = ifelse(j == 1, V[i, 2, k], δyᵃᶜᵃ(i, j, k, grid, V))

# Derivative Operators
@inline ∂xᶠᵃᵃ_c(i, j, k, grid, bc, η) = δxᶠᵃᵃ_c(i, j, k, grid, bc, η) / Δxᶠᶜᶜ(i, j, k, grid)
@inline ∂yᵃᶠᵃ_c(i, j, k, grid, bc, η) = δyᵃᶠᵃ_c(i, j, k, grid, bc, η) / Δyᶜᶠᶜ(i, j, k, grid)

# Derivative Operators
@inline ∂xᶜᵃᵃ_U(i, j, k, grid, bc, U) = δxᶜᵃᵃ_U(i, j, k, grid, bc, U) / Δxᶜᶜᶜ(i, j, k, grid)
@inline ∂yᵃᶜᵃ_V(i, j, k, grid, bc, V) = δyᵃᶜᵃ_V(i, j, k, grid, bc, V) / Δyᶜᶜᶜ(i, j, k, grid)
                            
const c = Center()
const f = Face()

@inline δx_f(ℓy, ℓz, i, j, k, ibg::IBG, bc, δx, args...) = 
    ifelse(immersed_inactive_node(i,   j, k, ibg, c, ℓy, ℓz) |
           immersed_inactive_node(i-1, j, k, ibg, c, ℓy, ℓz),
           zero(ibg),
           δx(i, j, k, ibg.underlying_grid, bc, args...))

@inline δx_c(ℓy, ℓz, i, j, k, ibg::IBG, bc, δx, args...) = 
    ifelse(immersed_inactive_node(i,   j, k, ibg, f, ℓy, ℓz) |
           immersed_inactive_node(i+1, j, k, ibg, f, ℓy, ℓz),
           zero(ibg),
           δx(i, j, k, ibg.underlying_grid, bc, args...))

@inline δy_f(ℓx, ℓz, i, j, k, ibg::IBG, bc, δy, args...) = 
    ifelse(immersed_inactive_node(i, j,   k, ibg, ℓx, c, ℓz) |
           immersed_inactive_node(i, j-1, k, ibg, ℓx, c, ℓz),
           zero(ibg),
           δy(i, j, k, ibg.underlying_grid, bc, args...))

@inline δy_c(ℓx, ℓz, i, j, k, ibg::IBG, bc, δy, args...) = 
    ifelse(immersed_inactive_node(i, j,   k, ibg, ℓx, f, ℓz) |
           immersed_inactive_node(i, j+1, k, ibg, ℓx, f, ℓz),
           zero(ibg),
           δy(i, j, k, ibg.underlying_grid, bc, args...))

@inline δx_f(ℓy, ℓz, i, j, k, ibg::IBG, bc::NoSlip, δx, c) = 
    @inbounds ifelse(immersed_inactive_node(i,   j, k, ibg, c, ℓy, ℓz),
                     2c[i-1, j, k],
              ifelse(immersed_inactive_node(i-1, j, k, ibg, c, ℓy, ℓz),
                     2c[i, j, k],
                     δx(i, j, k, ibg.underlying_grid, bc, c)))

@inline δx_c(ℓy, ℓz, i, j, k, ibg::IBG, bc::NoSlip, δx, U) = 
    @inbounds ifelse(immersed_inactive_node(i,   j, k, ibg, f, ℓy, ℓz),
                     2U[i+1, j, k],
              ifelse(immersed_inactive_node(i+1, j, k, ibg, f, ℓy, ℓz),
                     2U[i, j, k],
                     δx(i, j, k, ibg.underlying_grid, bc, U)))

@inline δy_f(ℓx, ℓz, i, j, k, ibg::IBG, bc::NoSlip, δy, c) = 
    @inbounds ifelse(immersed_inactive_node(i, j,   k, ibg, ℓx, c, ℓz),
                     2c[i, j-1, k],
              ifelse(immersed_inactive_node(i, j-1, k, ibg, ℓx, c, ℓz),
                     2c[i, j, k],
                     δy(i, j, k, ibg.underlying_grid, bc, c)))

@inline δy_c(ℓx, ℓz, i, j, k, ibg::IBG, bc::NoSlip, δy, V) = 
    @inbounds ifelse(immersed_inactive_node(i, j,   k, ibg, ℓx, f, ℓz),
                     2V[i, j+1, k],
              ifelse(immersed_inactive_node(i, j+1, k, ibg, ℓx, f, ℓz),
                     2V[i, j, k],
                     δy(i, j, k, ibg.underlying_grid, bc, V)))

@inline ∂xᶠᶠᶜ_c(i, j, k, grid, args...) = ∂xᶠᵃᵃ_c(i, j, k, grid, args...)
@inline ∂yᶠᶠᶜ_c(i, j, k, grid, args...) = ∂yᵃᶠᵃ_c(i, j, k, grid, args...)
@inline ∂xᶜᶠᶜ_U(i, j, k, grid, args...) = ∂xᶜᵃᵃ_U(i, j, k, grid, args...)
@inline ∂yᶠᶜᶜ_V(i, j, k, grid, args...) = ∂yᵃᶜᵃ_V(i, j, k, grid, args...)

@inline ∂xᶠᶠᶜ_c(i, j, k, grid::IBG, bc, c) = δx_f(f, c, i, j, k, grid, bc, ∂xᶠᵃᵃ_c, c)
@inline ∂yᶠᶠᶜ_c(i, j, k, grid::IBG, bc, c) = δy_f(f, c, i, j, k, grid, bc, ∂yᵃᶠᵃ_c, c)
@inline ∂xᶜᶠᶜ_U(i, j, k, grid::IBG, bc, U) = δx_c(f, c, i, j, k, grid, bc, ∂xᶜᵃᵃ_U, U)
@inline ∂yᶠᶜᶜ_V(i, j, k, grid::IBG, bc, V) = δy_c(f, c, i, j, k, grid, bc, ∂yᵃᶜᵃ_V, V)

@inline ∂xᶠᶜᶜ_c(i, j, k, grid, args...) = ∂xᶠᵃᵃ_c(i, j, k, grid, args...)
@inline ∂yᶜᶠᶜ_c(i, j, k, grid, args...) = ∂yᵃᶠᵃ_c(i, j, k, grid, args...)
@inline ∂xᶜᶜᶜ_U(i, j, k, grid, args...) = ∂xᶜᵃᵃ_U(i, j, k, grid, args...)
@inline ∂yᶜᶜᶜ_V(i, j, k, grid, args...) = ∂yᵃᶜᵃ_V(i, j, k, grid, args...)

@inline ∂xᶠᶜᶜ_c(i, j, k, grid::IBG, bc, c) = δx_f(c, c, i, j, k, grid, bc, ∂xᶠᵃᵃ_c, c)
@inline ∂yᶜᶠᶜ_c(i, j, k, grid::IBG, bc, c) = δy_f(c, c, i, j, k, grid, bc, ∂yᵃᶠᵃ_c, c)
@inline ∂xᶜᶜᶜ_U(i, j, k, grid::IBG, bc, U) = δx_c(c, c, i, j, k, grid, bc, ∂xᶜᵃᵃ_U, U)
@inline ∂yᶜᶜᶜ_V(i, j, k, grid::IBG, bc, V) = δy_c(c, c, i, j, k, grid, bc, ∂yᵃᶜᵃ_V, V)