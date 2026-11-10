#####
##### Boundary contributions to the column energy solve. External heat fluxes are upward-positive: a direct flux enters
##### the top cell as -Qᵘ and the bottom cell as +Qᵇ. Dirichlet boundaries add an implicit conductance and a conductive flux.
#####

@inline function column_external_flux_temperature_slope(i, j, grid, ext, Ts, clock, model_fields)
    δT = cbrt(eps(eltype(grid))) * max(1, abs(Ts))
    Qplus = getflux(ext, i, j, grid, Ts + δT, clock, model_fields)
    Qminus = getflux(ext, i, j, grid, Ts - δT, clock, model_fields)
    return (Qplus - Qminus) / (2δT)
end

@inline function column_external_flux_temperature_slope(i, j, grid, ext::LinearizedSurfaceFlux, Ts, clock, model_fields)
    return getflux(ext.temperature_derivative, i, j, grid, Ts, clock, model_fields)
end

@inline function column_bottom_boundary_energy_factor(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    λ = column_external_flux_temperature_slope(i, j, grid, ext, T, clock, model_fields)
    return -Δt * λ * fields.temperature_energy_derivative[i, j, k] / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_top_boundary_energy_factor(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    λ = column_external_flux_temperature_slope(i, j, grid, ext, T, clock, model_fields)
    return Δt * λ * fields.temperature_energy_derivative[i, j, k] / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_bottom_boundary_energy_factor(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_top_boundary_energy_factor(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_bottom_boundary_energy_flux(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    λ = column_external_flux_temperature_slope(i, j, grid, ext, T, clock, model_fields)
    return - getflux(ext, i, j, grid, T, clock, model_fields) + λ * fields.temperature_energy_derivative[i, j, k] * fields.internal_energy[i, j, k]
end

@inline function column_top_boundary_energy_flux(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    λ = column_external_flux_temperature_slope(i, j, grid, ext, T, clock, model_fields)
    return - getflux(ext, i, j, grid, T, clock, model_fields) + λ * fields.temperature_energy_derivative[i, j, k] * fields.internal_energy[i, j, k]
end

# A Dirichlet base couples to the interior by conduction only; its external flux enters the basal Stefan balance.
@inline function column_bottom_boundary_energy_flux(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(i, j, grid, bc, relation, model_fields)
    return conductance * (T - Tᵇ - TE * E)
end

@inline function column_top_boundary_energy_flux(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(i, j, grid, bc, relation, model_fields)
    return conductance * (Tᵇ - T + TE * E) - getflux(ext, i, j, grid, Tᵇ, clock, model_fields)
end

@inline column_surface_stefan_residual_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(eltype(grid))

@inline function column_surface_stefan_residual_flux(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Δz = Δzᶜᶜᶜ(i, j, k, grid)
    # the cap is evaluated against the start-of-step enthalpy, the same one the surface solve used
    E  = @inbounds auxiliary.surface_start_energy[i, j, 1]
    S  = @inbounds fields.bulk_salinity[i, j, k]
    Tₛ = @inbounds auxiliary.surface_temperature[i, j, 1]
    requested = - getflux(ext, i, j, grid, Tₛ, clock, model_fields)
    available = column_energy_to_complete_melt(relation, E, S, Δz) / Δt
    applied = min(requested, max(available, zero(available)))
    return applied - requested
end

