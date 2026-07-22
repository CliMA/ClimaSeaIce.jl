using Oceananigans.Advection: FluxFormAdvection,
                              _advective_tracer_flux_x,
                              _advective_tracer_flux_y,
                              conditional_flux_fcc,
                              conditional_flux_cfc

# To obtain better numerical properties (conservation of `ℵ·h` in convergent flow
# in the presence of WENO over/undershoots), ice volume per area `V = ℵ·h` is
# transported instead of the bare thickness `h`. The tendency
#
#     ∂V/∂t = -∇·(U·V)
#
# is computed below by `div_Uℵh`. The recovery `h = V/ℵ` is then performed in
# `_dynamic_step_tracers!` with a small-ℵ guard, following SI³
# (Vancoppenolle et al. 2020, `icedyn_adv_pra.F90`).

_advective_thickness_flux_x(i, j, k, grid, scheme, U, ℵ, h) = advective_thickness_flux_x(i, j, k, grid, scheme, U, ℵ, h)
_advective_thickness_flux_y(i, j, k, grid, scheme, U, ℵ, h) = advective_thickness_flux_y(i, j, k, grid, scheme, U, ℵ, h)

_advective_thickness_flux_x(i, j, k, ibg::ImmersedBoundaryGrid, scheme, U, ℵ, h) =
    conditional_flux_fcc(i, j, k, ibg, zero(ibg), advective_thickness_flux_x(i, j, k, ibg, scheme, U, ℵ, h))

_advective_thickness_flux_y(i, j, k, ibg::ImmersedBoundaryGrid, scheme, U, ℵ, h) =
    conditional_flux_cfc(i, j, k, ibg, zero(ibg), advective_thickness_flux_y(i, j, k, ibg, scheme, U, ℵ, h))

@inline function advective_thickness_flux_x(i, j, k, grid, advection, U, ℵ, h)
    ϕℵ = advective_tracer_flux_x(i, j, k, grid, advection, U, ℵ) / Axᶠᶜᶜ(i, j, k, grid)
    Uϕℵh = ϕℵ * advective_tracer_flux_x(i, j, k, grid, advection, U, h)
    @inbounds ϕℵh = ifelse(U[i, j, k] == 0, zero(grid), Uϕℵh / U[i, j, k])
    return ϕℵh
end

@inline function advective_thickness_flux_y(i, j, k, grid, advection, V, ℵ, h)
    ϕℵ = advective_tracer_flux_y(i, j, k, grid, advection, V, ℵ) / Ayᶜᶠᶜ(i, j, k, grid)
    Vϕℵh = ϕℵ * advective_tracer_flux_y(i, j, k, grid, advection, V, h)
    @inbounds ϕℵh = ifelse(V[i, j, k] == 0, zero(grid), Vϕℵh / V[i, j, k])
    return ϕℵh
end

@inline div_Uℵh(i, j, k, grid, ::Nothing, U, ℵ, h) = zero(grid)

# Volume-per-area tendency: returns ∇·(U·ℵ·h).
@inline function div_Uℵh(i, j, k, grid, advection, U, ℵ, h)
    return 1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_thickness_flux_x, advection, U.u, ℵ, h) +
                                      δyᵃᶜᵃ(i, j, k, grid, _advective_thickness_flux_y, advection, U.v, ℵ, h))
end

# Volume-per-area tendency: returns ∇·(U·ℵ·h).
@inline function div_Uℵh(i, j, k, grid, advection::FluxFormAdvection, U, ℵ, h)
    return 1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_thickness_flux_x, advection.x, U.u, ℵ, h) +
                                      δyᵃᶜᵃ(i, j, k, grid, _advective_thickness_flux_y, advection.y, U.v, ℵ, h))
end

@inline horizontal_div_Uc(i, j, k, grid, ::Nothing, U, c) = zero(grid)
@inline horizontal_div_Uc(i, j, k, grid, advection, U, c) =
    1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_tracer_flux_x, advection, U.u, c) +
                               δyᵃᶜᵃ(i, j, k, grid, _advective_tracer_flux_y, advection, U.v, c))

@inline horizontal_div_Uc(i, j, k, grid, advection::FluxFormAdvection, U, c) =
    1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_tracer_flux_x, advection.x, U.u, c) +
                               δyᵃᶜᵃ(i, j, k, grid, _advective_tracer_flux_y, advection.y, U.v, c))
