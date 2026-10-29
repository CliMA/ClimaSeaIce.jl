#####
##### Boundary contributions to the column energy solve. External heat fluxes are upward-positive: a direct flux enters
##### the top cell as -Qᵘ and the bottom cell as +Qᵇ. Dirichlet boundaries add an implicit conductance and a conductive flux.
#####

@inline column_bottom_boundary_energy_factor(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)
@inline column_top_boundary_energy_factor(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)

# On a cold surface, Qext(Ts) = C(Ttop - Ts). Differentiating this
# massless-surface balance gives dQext/dTtop = C Q′/(C + Q′). This Robin
# conductance belongs in the same tridiagonal solve as interior conduction.
@inline function column_external_flux_temperature_slope(ext, i, j, grid, Ts, clock, model_fields)
    δT = oftype(Ts, 1e-3)
    Qplus = getflux(ext, i, j, grid, Ts + δT, clock, model_fields)
    Qminus = getflux(ext, i, j, grid, Ts - δT, clock, model_fields)
    return max(zero(Ts), (Qplus - Qminus) / (2δT))
end

@inline column_external_flux_temperature_slope(ext::Number, i, j, grid, Ts, clock, model_fields) = zero(Ts)
@inline column_external_flux_temperature_slope(ext::AbstractArray, i, j, grid, Ts, clock, model_fields) = zero(Ts)

@inline function column_external_flux_temperature_slope(ext::LinearizedSurfaceFlux, i, j, grid, Ts, clock, model_fields)
    return max(zero(Ts), getflux(ext.temperature_derivative, i, j, grid, Ts, clock, model_fields))
end

@inline function column_cold_surface_conductance(ext, i, j, k, grid, auxiliary, clock, model_fields)
    C = column_boundary_temperature_conductance(i, j, size(grid, 3)+1, k, grid, auxiliary)
    Ts = @inbounds auxiliary.surface_temperature[i, j, 1]
    Q′ = column_external_flux_temperature_slope(ext, i, j, grid, Ts, clock, model_fields)
    return C * Q′ / (C + Q′)
end

@inline function column_surface_is_melting(i, j, k, auxiliary, fields, relation)
    @inbounds begin
        Tₛ = auxiliary.surface_temperature[i, j, 1]
        Tₘ = melting_temperature(relation.phase_transitions.liquidus, fields.bulk_salinity[i, j, k])
    end
    return Tₛ >= Tₘ
end

@inline function column_bottom_boundary_energy_factor(bc::ColumnDirichletBoundary, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_top_boundary_energy_factor(bc::ColumnDirichletBoundary, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_top_boundary_energy_factor(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = ifelse(column_surface_is_melting(i, j, Nz, auxiliary, fields, relation),
                         column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary),
                         column_cold_surface_conductance(ext, i, j, k, grid, auxiliary, clock, model_fields))
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_bottom_boundary_energy_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    return - getflux(ext, i, j, grid, T, clock, model_fields)
end

@inline function column_top_boundary_energy_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    return - getflux(ext, i, j, grid, T, clock, model_fields)
end

# A Dirichlet base couples to the interior by conduction only; its external flux enters the basal Stefan balance.
@inline function column_bottom_boundary_energy_flux(bc::ColumnDirichletBoundary, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(bc, i, j, grid, relation)
    return conductance * (T - Tᵇ - TE * E)
end

@inline function column_top_boundary_energy_flux(bc::ColumnDirichletBoundary, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(bc, i, j, grid, relation)
    return conductance * (Tᵇ - T + TE * E) - getflux(ext, i, j, grid, Tᵇ, clock, model_fields)
end

# At a melting surface, only the conductive flux enters the top cell. The
# atmospheric energy in excess of conduction ablates ice at the surface.
@inline function column_top_boundary_energy_flux(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Tₛ = @inbounds auxiliary.surface_temperature[i, j, 1]
    cold_conductance = column_cold_surface_conductance(ext, i, j, k, grid, auxiliary, clock, model_fields)
    conductance = column_boundary_temperature_conductance(i, j, size(grid, 3)+1, k, grid, auxiliary)
    @inbounds begin
        T = fields.temperature[i, j, k]
        E = fields.internal_energy[i, j, k]
        TE = fields.temperature_energy_derivative[i, j, k]
    end
    cold_flux = -getflux(ext, i, j, grid, Tₛ, clock, model_fields) + cold_conductance * TE * E
    melting_flux = conductance * (Tₛ - T + TE * E)
    return ifelse(column_surface_is_melting(i, j, k, auxiliary, fields, relation), melting_flux, cold_flux)
end

@inline column_surface_stefan_residual_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(eltype(grid))

@inline function column_surface_stefan_residual_flux(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Tₛ = @inbounds auxiliary.surface_temperature[i, j, 1]
    T = @inbounds fields.temperature[i, j, k]
    conductance = column_boundary_temperature_conductance(i, j, k+1, k, grid, auxiliary)
    Qᶜ = conductance * (T - Tₛ)
    Qᵃ = getflux(ext, i, j, grid, Tₛ, clock, model_fields)
    return ifelse(column_surface_is_melting(i, j, k, auxiliary, fields, relation),
                  min(Qᵃ - Qᶜ, zero(Qᵃ)), zero(Qᵃ))
end

