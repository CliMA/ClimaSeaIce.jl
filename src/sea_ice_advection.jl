using Oceananigans.Operators
using Oceananigans.ImmersedBoundaries: IBG
using Oceananigans.Advection: FluxFormAdvection, 
                              _advective_tracer_flux_x, 
                              _advective_tracer_flux_y,
                              conditional_flux_fcc, 
                              conditional_flux_cfc,
                              UpwindScheme, 
                              CenteredScheme

using Oceananigans.Advection: _biased_interpolate_xᶠᵃᵃ, _biased_interpolate_yᵃᶠᵃ, bias

@inline _advective_sea_ice_tracer_flux_x(i, j, k, grid, scheme, u, c) = advective_sea_ice_tracer_flux_x(i, j, k, grid, scheme, u, c)
@inline _advective_sea_ice_tracer_flux_y(i, j, k, grid, scheme, v, c) = advective_sea_ice_tracer_flux_y(i, j, k, grid, scheme, v, c)

@inline _advective_sea_ice_tracer_flux_x(i, j, k, ibg::IBG, scheme, u, c) = conditional_flux_fcc(i, j, k, ibg, zero(ibg), advective_sea_ice_tracer_flux_x(i, j, k, ibg, scheme, u, c))
@inline _advective_sea_ice_tracer_flux_y(i, j, k, ibg::IBG, scheme, v, c) = conditional_flux_cfc(i, j, k, ibg, zero(ibg), advective_sea_ice_tracer_flux_y(i, j, k, ibg, scheme, v, c))

@inline function advective_sea_ice_tracer_flux_x(i, j, k, grid, scheme::UpwindScheme, u, c)
    ũ  = ℑyᵃᶜᵃ(i, j, k, grid, Ax_qᶠᶠᶜ, u)
    cᴿ = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, bias(ũ), c)
    return ũ * cᴿ
end

@inline function advective_sea_ice_tracer_flux_y(i, j, k, grid, scheme::UpwindScheme, v, c)
    ṽ  = ℑxᶜᵃᵃ(i, j, k, grid, Ay_qᶠᶠᶜ, v)
    cᴿ = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, bias(ṽ), c)
    return ṽ * cᴿ
end

@inline function advective_sea_ice_tracer_flux_x(i, j, k, grid, scheme::CenteredScheme, u, c)
    ũ  = ℑyᵃᶜᵃ(i, j, k, grid, Ax_qᶠᶠᶜ, u)
    cᴿ = _symmetric_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, c)
    return ũ * cᴿ
end

@inline function advective_sea_ice_tracer_flux_y(i, j, k, grid, scheme::CenteredScheme, v, c)
    ṽ  = ℑxᶜᵃᵃ(i, j, k, grid, Ay_qᶠᶠᶜ, v)
    cᴿ = _symmetric_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, c)
    return ṽ * cᴿ
end

@inline horizontal_div_Uc(i, j, k, grid, ::Nothing, U, c) = zero(grid)

@inline horizontal_div_Uc(i, j, k, grid, scheme, U, c) = 
    1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_sea_ice_tracer_flux_x, scheme, U.u, c) +
                               δyᵃᶜᵃ(i, j, k, grid, _advective_sea_ice_tracer_flux_y, scheme, U.v, c))
                               
@inline horizontal_div_Uc(i, j, k, grid, scheme::FluxFormAdvection, U, c) = 
    1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_sea_ice_tracer_flux_x, scheme.x, U.u, c) +
                               δyᵃᶜᵃ(i, j, k, grid, _advective_sea_ice_tracer_flux_y, scheme.y, U.v, c))
                               
