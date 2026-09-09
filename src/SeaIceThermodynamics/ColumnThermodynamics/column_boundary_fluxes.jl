#####
##### Boundary contributions to the column energy solve, keyed on (heat_boundary_condition, external_flux).
##### `external_heat_fluxes` use the same upward-positive convention as `SlabThermodynamics` (positive = heat
##### leaving the ice, i.e. into the air at the top and from the ocean into the base at the bottom), so the same
##### `model.external_heat_fluxes` drives a slab and a column identically. A direct flux therefore enters the top
##### cell as `-Qᵘ` and the bottom cell as `+Qᵇ`. Dirichlet boundaries (`PrescribedTemperature`,
##### `IceWaterThermalEquilibrium`) additionally add an implicit conductance factor (LHS) and a conductive flux.
#####

@inline column_bottom_boundary_energy_factor(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)
@inline column_top_boundary_energy_factor(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)

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

# A direct-flux boundary injects its paired upward-positive external flux: a `left_flux` of `-Qᵇ` adds ocean heat
# to the base, a `right_flux` of `-Qᵘ` removes heat to the air at the top.
@inline function column_bottom_boundary_energy_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    return - getflux(ext, i, j, grid, T, clock, model_fields)
end

@inline function column_top_boundary_energy_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    return - getflux(ext, i, j, grid, T, clock, model_fields)
end

# A Dirichlet base pins the interface temperature, so the interior couples to it by conduction only — the paired
# ocean/lake `external_heat_fluxes` enters the basal Stefan balance (`column_basal_stefan_flux`), not the interior
# solve, exactly as in SlabThermodynamics (and CICE: basal growth = (Q_cond − Q_ocean)/ℒ).
@inline function column_bottom_boundary_energy_flux(bc::ColumnDirichletBoundary,
                                                    ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
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

# A `MeltingConstrainedFluxBalance` surface injects the downward flux `-Qᵘ(Tₛ)` evaluated at the iteratively-solved
# surface temperature `Tₛ`, capped at the top-cell complete-melt energy; the excess becomes the surface Stefan
# residual that drives surface melt.
@inline function column_top_boundary_energy_flux(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Δz = Δzᶜᶜᶜ(i, j, k, grid)
    E  = @inbounds fields.internal_energy[i, j, k]
    S  = @inbounds fields.bulk_salinity[i, j, k]
    Tₛ = @inbounds auxiliary.surface_temperature[i, j, 1]
    requested = - getflux(ext, i, j, grid, Tₛ, clock, model_fields)
    available = column_energy_to_complete_melt(relation, E, S, Δz) / Δt
    return min(requested, max(available, zero(available)))
end

# The surface Stefan residual is the part of the requested top flux not absorbed by the capped cell warming.
@inline column_surface_stefan_residual_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(eltype(grid))

@inline function column_surface_stefan_residual_flux(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Δz = Δzᶜᶜᶜ(i, j, k, grid)
    # Use the start-of-step enthalpy that the surface solve capped against, not the post-solve enthalpy: the
    # capped warming already drove the cell toward complete melt, so reading the current E would compute a near-zero
    # `available` and ablate the surface against the full requested flux (double-counting the absorbed energy).
    E  = @inbounds auxiliary.surface_start_energy[i, j, 1]
    S  = @inbounds fields.bulk_salinity[i, j, k]
    Tₛ = @inbounds auxiliary.surface_temperature[i, j, 1]
    requested = - getflux(ext, i, j, grid, Tₛ, clock, model_fields)
    available = column_energy_to_complete_melt(relation, E, S, Δz) / Δt
    applied = min(requested, max(available, zero(available)))
    return applied - requested
end

