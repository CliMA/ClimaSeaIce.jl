#####
##### Boundary contributions to the column energy solve. External heat fluxes are upward-positive: a direct flux enters
##### the top cell as -Qᵘ and the bottom cell as +Qᵇ. Dirichlet boundaries add an implicit conductance and a conductive flux.
#####

@inline column_bottom_boundary_energy_factor(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)
@inline column_top_boundary_energy_factor(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)

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
    column_surface_is_melting(i, j, Nz, auxiliary, fields, relation) || return zero(Δt)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
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
    column_surface_is_melting(i, j, k, auxiliary, fields, relation) ||
        return -getflux(ext, i, j, grid, Tₛ, clock, model_fields)
    conductance = column_boundary_temperature_conductance(i, j, size(grid, 3)+1, k, grid, auxiliary)
    @inbounds begin
        T = fields.temperature[i, j, k]
        E = fields.internal_energy[i, j, k]
        TE = fields.temperature_energy_derivative[i, j, k]
    end
    return conductance * (Tₛ - T + TE * E)
end

@inline column_surface_stefan_residual_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(eltype(grid))

@inline function column_surface_stefan_residual_flux(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    column_surface_is_melting(i, j, k, auxiliary, fields, relation) || return zero(Δt)
    Tₛ = @inbounds auxiliary.surface_temperature[i, j, 1]
    T = @inbounds fields.temperature[i, j, k]
    conductance = column_boundary_temperature_conductance(i, j, k+1, k, grid, auxiliary)
    Qᶜ = conductance * (T - Tₛ)
    Qᵃ = getflux(ext, i, j, grid, Tₛ, clock, model_fields)
    return min(Qᵃ - Qᶜ, zero(Qᵃ))
end

