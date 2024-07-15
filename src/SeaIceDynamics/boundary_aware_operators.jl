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
@inline ı(i, j, k, grid, ϕ)                    = ϕ[i, j, k]

# Defining Interpolation operators for the immersed boundaries
@inline conditional_ℑx_f(LY, LZ, i, j, k, aug::AUG, ℑx, args...) = ℑx(i, j, k, aug, args...)
@inline conditional_ℑx_c(LY, LZ, i, j, k, aug::AUG, ℑx, args...) = ℑx(i, j, k, aug, args...)
@inline conditional_ℑy_f(LX, LZ, i, j, k, aug::AUG, ℑy, args...) = ℑy(i, j, k, aug, args...)
@inline conditional_ℑy_c(LX, LZ, i, j, k, aug::AUG, ℑy, args...) = ℑy(i, j, k, aug, args...)

@inline conditional_ℑx_f(LY, LZ, i, j, k, ibg::IBG, ℑx, args...) = ifelse(inactive_node(i, j, k, ibg, c, LY, LZ), ı(i-1, j, k, ibg, args...), ifelse(inactive_node(i-1, j, k, ibg, c, LY, LZ), ı(i, j, k, ibg, args...), ℑx(i, j, k, ibg.underlying_grid, args...)))
@inline conditional_ℑx_c(LY, LZ, i, j, k, ibg::IBG, ℑx, args...) = ifelse(inactive_node(i, j, k, ibg, f, LY, LZ), ı(i+1, j, k, ibg, args...), ifelse(inactive_node(i+1, j, k, ibg, f, LY, LZ), ı(i, j, k, ibg, args...), ℑx(i, j, k, ibg.underlying_grid, args...)))
@inline conditional_ℑy_f(LX, LZ, i, j, k, ibg::IBG, ℑy, args...) = ifelse(inactive_node(i, j, k, ibg, LX, c, LZ), ı(i, j-1, k, ibg, args...), ifelse(inactive_node(i, j-1, k, ibg, LX, c, LZ), ı(i, j, k, ibg, args...), ℑy(i, j, k, ibg.underlying_grid, args...)))
@inline conditional_ℑy_c(LX, LZ, i, j, k, ibg::IBG, ℑy, args...) = ifelse(inactive_node(i, j, k, ibg, LX, f, LZ), ı(i, j+1, k, ibg, args...), ifelse(inactive_node(i, j+1, k, ibg, LX, f, LZ), ı(i, j, k, ibg, args...), ℑy(i, j, k, ibg.underlying_grid, args...)))

@inline ℑxᴮᶜᶜᶜ(i, j, k, ibg, args...) = conditional_ℑx_c(c, c, i, j, k, ibg, ℑxᶜᵃᵃ, args...)
@inline ℑxᴮᶠᶜᶜ(i, j, k, ibg, args...) = conditional_ℑx_f(c, c, i, j, k, ibg, ℑxᶠᵃᵃ, args...)
@inline ℑyᴮᶜᶜᶜ(i, j, k, ibg, args...) = conditional_ℑy_c(c, c, i, j, k, ibg, ℑyᵃᶜᵃ, args...)
@inline ℑyᴮᶜᶠᶜ(i, j, k, ibg, args...) = conditional_ℑy_f(c, c, i, j, k, ibg, ℑyᵃᶠᵃ, args...)

@inline ℑxᴮᶜᶠᶜ(i, j, k, ibg, args...) = conditional_ℑx_c(f, c, i, j, k, ibg, ℑxᶜᵃᵃ, args...)
@inline ℑxᴮᶠᶠᶜ(i, j, k, ibg, args...) = conditional_ℑx_f(f, c, i, j, k, ibg, ℑxᶠᵃᵃ, args...)
@inline ℑyᴮᶠᶜᶜ(i, j, k, ibg, args...) = conditional_ℑy_c(f, c, i, j, k, ibg, ℑyᵃᶜᵃ, args...)
@inline ℑyᴮᶠᶠᶜ(i, j, k, ibg, args...) = conditional_ℑy_f(f, c, i, j, k, ibg, ℑyᵃᶠᵃ, args...)

@inline ℑxyᴮᶜᶜᶜ(i, j, k, grid, args...) = ℑxᴮᶜᶜᶜ(i, j, k, grid, ℑyᴮᶜᶜᶜ, args...)
@inline ℑxyᴮᶠᶜᶜ(i, j, k, grid, args...) = ℑxᴮᶠᶜᶜ(i, j, k, grid, ℑyᴮᶠᶜᶜ, args...)
@inline ℑxyᴮᶜᶠᶜ(i, j, k, grid, args...) = ℑxᴮᶜᶠᶜ(i, j, k, grid, ℑyᴮᶜᶠᶜ, args...)
@inline ℑxyᴮᶠᶠᶜ(i, j, k, grid, args...) = ℑxᴮᶠᶠᶜ(i, j, k, grid, ℑyᴮᶠᶠᶜ, args...)