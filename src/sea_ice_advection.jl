using Oceananigans.Operators
using Oceananigans.ImmersedBoundaries
using Oceananigans.Advection: FluxFormAdvection, 
                              _advective_tracer_flux_x, 
                              _advective_tracer_flux_y,
                              conditional_flux_fcc, 
                              conditional_flux_cfc

# To obtain better numerical properties, the ice thickness is advected together 
# with the concentration, i.e.:
# 
# A = ℵ⁻¹ ∇ ⋅ (uℵh) 
#
# instead of the classical 
# 
# A = ∇ ⋅ (uh)

using Oceananigans.Advection: _biased_interpolate_xᶠᵃᵃ, _biased_interpolate_yᵃᶠᵃ, bias

@inline function advective_new_tracer_flux_x(i, j, k, grid, advection, u, c)
    ũ  = ℑyᵃᶜᵃ(i, j, k, grid, u)
    cᴿ = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, advection, bias(ũ), c)
    return Axᶠᶜᶜ(i, j, k, grid) * ũ * cᴿ
end

@inline function advective_new_tracer_flux_y(i, j, k, grid, advection, v, c)
    ṽ  = ℑxᶜᵃᵃ(i, j, k, grid, v)
    cᴿ = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, advection, bias(ṽ), c)
    return Ayᶜᶠᶜ(i, j, k, grid) * ṽ * cᴿ
end

@inline horizontal_div_Uc(i, j, k, grid, ::Nothing, U, c) = zero(grid)
@inline horizontal_div_Uc(i, j, k, grid, advection, U, c) = 
    1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, advective_new_tracer_flux_x, advection, U.u, c) +
                               δyᵃᶜᵃ(i, j, k, grid, advective_new_tracer_flux_y, advection, U.v, c))
                               
@inline horizontal_div_Uc(i, j, k, grid, advection::FluxFormAdvection, U, c) = 
    1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, advective_new_tracer_flux_x, advection.x, U.u, c) +
                               δyᵃᶜᵃ(i, j, k, grid, advective_new_tracer_flux_y, advection.y, U.v, c))
                               
