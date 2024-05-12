using Oceananigans.Operators
using Oceananigans.Grids: peripheral_node

# Identity function for the interpolation operators
@inline ı(i, j, k, grid, f::Function, args...) = f(i, j, k, grid, args...)
@inline ı(i, j, k, grid, ϕ)                    = @inbounds ϕ[i, j, k]

# Defining Interpolation operators for the immersed boundaries
@inline conditional_ℑx_f(LY, LZ, i, j, k, grid, t, args...) = ifelse(peripheral_node(i,   j, k, grid, Center(), LY, LZ), ı(i-1, j, k, grid, t, args...), 
                                                              ifelse(peripheral_node(i-1, j, k, grid, Center(), LY, LZ), ı(i,   j, k, grid, t, args...), 
                                                              ℑxᶠᵃᵃ(i, j, k, grid, t, args...)))

@inline conditional_ℑy_f(LX, LZ, i, j, k, grid, t, args...) = ifelse(peripheral_node(i, j,   k, grid, LX, Center(), LZ), ı(i, j-1, k, grid, t, args...), 
                                                              ifelse(peripheral_node(i, j-1, k, grid, LX, Center(), LZ), ı(i, j,   k, grid, t, args...), 
                                                              ℑyᵃᶠᵃ(i, j, k, grid, t, args...)))

@inline conditional_ℑx_c(LY, LZ, i, j, k, grid, t, args...) = ifelse(peripheral_node(i,   j, k, grid, Face(), LY, LZ), ı(i+1, j, k, grid, t, args...), 
                                                              ifelse(peripheral_node(i+1, j, k, grid, Face(), LY, LZ), ı(i,   j, k, grid, t, args...), 
                                                              ℑxᶜᵃᵃ(i, j, k, grid, t, args...)))

@inline conditional_ℑy_c(LX, LZ, i, j, k, grid, t, args...) = ifelse(peripheral_node(i, j,   k, grid, LX, Face(), LZ), ı(i, j+1, k, grid, t, args...), 
                                                              ifelse(peripheral_node(i, j+1, k, grid, LX, Face(), LZ), ı(i, j,   k, grid, t, args...),
                                                               ℑyᵃᶜᵃ(i, j, k, grid, t, args...)))

ℑxᶠᶜᶜ(i, j, k, grid, t, args...) = conditional_ℑx_f(Center(), Center(), i, j, k, grid, t, args...)
ℑyᶜᶠᶜ(i, j, k, grid, t, args...) = conditional_ℑy_f(Center(), Center(), i, j, k, grid, t, args...)

ℑxᶠᶠᶜ(i, j, k, grid, t, args...) = conditional_ℑx_f(Face(), Center(), i, j, k, grid, t, args...)
ℑyᶠᶠᶜ(i, j, k, grid, t, args...) = conditional_ℑy_f(Face(), Center(), i, j, k, grid, t, args...)

ℑxᶜᶜᶜ(i, j, k, grid, t, args...) = conditional_ℑx_c(Center(), Center(), i, j, k, grid, t, args...)
ℑyᶜᶜᶜ(i, j, k, grid, t, args...) = conditional_ℑy_c(Center(), Center(), i, j, k, grid, t, args...)

# The only four we need?
ℑxyᶠᶠᶜ(i, j, k, grid, t, args...) = ℑxᶠᶠᶜ(i, j, k, grid, ℑyᶜᶠᶜ, t, args...)
ℑxyᶜᶜᶜ(i, j, k, grid, t, args...) = ℑyᶠᶠᶜ(i, j, k, grid, ℑxᶠᶜᶜ, t, args...)
ℑxyᶠᶜᶜ(i, j, k, grid, t, args...) = ℑxᶠᶜᶜ(i, j, k, grid, ℑyᶜᶜᶜ, t, args...)
ℑxyᶜᶠᶜ(i, j, k, grid, t, args...) = ℑyᶜᶠᶜ(i, j, k, grid, ℑxᶜᶜᶜ, t, args...)