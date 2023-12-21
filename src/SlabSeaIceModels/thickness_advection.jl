using Oceananigans.Operators
using Oceananigans.Advection: _advective_tracer_flux_x, _advective_tracer_flux_y

@inline function _advective_thickness_flux_x(i, j, k, grid, advection, U, ℵ, h)
    ϕℵ = _advective_tracer_flux_x(i, j, k, grid, advection, U, ℵ) / Axᶠᶜᶜ(i, j, k, grid)
    @inbounds ϕh = ϕℵ * _advective_tracer_flux_x(i, j, k, grid, advection, U, h) / U[i, j, k] 
    return ϕh
end

@inline function _advective_thickness_flux_y(i, j, k, grid, advection, V, ℵ, h)
    ϕℵ = _advective_tracer_flux_x(i, j, k, grid, advection, V, ℵ) / Ayᶜᶠᶜ(i, j, k, grid)
    @inbounds ϕh = ϕℵ * _advective_tracer_flux_y(i, j, k, grid, advection, V, h) / V[i, j, k] 
    return ϕh
end

# For thickness, we compute [ℵ⁻¹ ∇ ⋅ (uℵh)]
@inline function div_Uℵh(i, j, grid, advection, U, ℵ, h)
    div_Uℵh = 1 / Vᶜᶜᶜ(i, j, 1, grid) * (δxᶜᵃᵃ(i, j, 1, grid, _advective_thickness_flux_x, advection, U.u, ℵ, h) +
                                         δyᵃᶜᵃ(i, j, 1, grid, _advective_thickness_flux_x, advection, U.v, ℵ, h))

    @inbounds ℵ⁻¹ = ifelse(ℵ[i, j, 1] > 1e-10, 1 / ℵ[i, j, 1], 0)

    return ℵ⁻¹ * div_Uℵh
end
