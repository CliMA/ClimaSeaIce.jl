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

_advective_thickness_flux_x(args...) = advective_thickness_flux_x(args...)
_advective_thickness_flux_y(args...) = advective_thickness_flux_y(args...)

_advective_thickness_flux_x(i, j, k, ibg::ImmersedBoundaryGrid, args...) = 
    conditional_flux_fcc(i, j, k, ibg, zero(ibg), advective_thickness_flux_x(i, j, k, ibg, args...))

_advective_thickness_flux_y(i, j, k, ibg::ImmersedBoundaryGrid, args...) = 
    conditional_flux_cfc(i, j, k, ibg, zero(ibg), advective_thickness_flux_y(i, j, k, ibg, args...))

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

@inline div_Uℵh(i, j, k, grid, ::Nothing, U, ℵ, h) = zero(grid)>>>>>>> main

# For thickness, we compute [ℵ⁻¹ ∇ ⋅ (uℵh)]
@inline function div_Uℵh(i, j, k, grid, advection, U, ℵ, h)
    ∇Uℵh = 1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_thickness_flux_x, advection, U.u, ℵ, h) +
                                      δyᵃᶜᵃ(i, j, k, grid, _advective_thickness_flux_y, advection, U.v, ℵ, h))

    @inbounds ℵ⁻¹ = ifelse(ℵ[i, j, k] != 0, 1 / ℵ[i, j, k], zero(grid))

    return ℵ⁻¹ * ∇Uℵh
end

# For thickness, we compute [ℵ⁻¹ ∇ ⋅ (uℵh)]
@inline function div_Uℵh(i, j, k, grid, advection::FluxFormAdvection, U, ℵ, h)
    ∇Uℵh = 1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_thickness_flux_x, advection.x, U.u, ℵ, h) +
                                      δyᵃᶜᵃ(i, j, k, grid, _advective_thickness_flux_y, advection.y, U.v, ℵ, h))

    @inbounds ℵ⁻¹ = ifelse(ℵ[i, j, k] != 0, 1 / ℵ[i, j, k], zero(grid))

    return ℵ⁻¹ * ∇Uℵh
end

@inline horizontal_div_Uc(i, j, k, grid, ::Nothing, U, c) = zero(grid)
@inline horizontal_div_Uc(i, j, k, grid, advection, U, c) = 
    1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_tracer_flux_x, advection, U.u, c) +
                               δyᵃᶜᵃ(i, j, k, grid, _advective_tracer_flux_y, advection, U.v, c))
                               
@inline horizontal_div_Uc(i, j, k, grid, advection::FluxFormAdvection, U, c) = 
    1 / Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_tracer_flux_x, advection.x, U.u, c) +
                               δyᵃᶜᵃ(i, j, k, grid, _advective_tracer_flux_y, advection.y, U.v, c))
                               
