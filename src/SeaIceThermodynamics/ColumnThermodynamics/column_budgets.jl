"""
    column_integrated_energy(thermodynamics)

Return the vertical integral of volumetric internal energy over every column
cell in `thermodynamics`.
"""
function column_integrated_energy(thermodynamics::ColumnEnergyThermodynamics)
    E = thermodynamics.fields.internal_energy
    grid = E.grid
    total_energy = zero(eltype(grid))

    for k in 1:size(grid, 3), j in 1:size(grid, 2), i in 1:size(grid, 1)
        total_energy += @inbounds E[i, j, k] * Δzᶜᶜᶜ(i, j, k, grid)
    end

    return total_energy
end

"""
    column_integrated_salinity(thermodynamics)

Return the vertical integral of bulk salinity over every column cell in
`thermodynamics`.
"""
function column_integrated_salinity(thermodynamics::ColumnEnergyThermodynamics)
    S = thermodynamics.fields.bulk_salinity
    grid = S.grid
    total_salinity = zero(eltype(grid))

    for k in 1:size(grid, 3), j in 1:size(grid, 2), i in 1:size(grid, 1)
        total_salinity += @inbounds S[i, j, k] * Δzᶜᶜᶜ(i, j, k, grid)
    end

    return total_salinity
end

function column_boundary_energy_flux_difference(thermodynamics::ColumnEnergyThermodynamics,
                                                external_heat_fluxes, clock, model_fields)
    grid = thermodynamics.fields.internal_energy.grid
    temperature = thermodynamics.fields.temperature
    Nz = size(grid, 3)
    difference = zero(eltype(grid))

    # Upward-positive fluxes: net energy gained by the column is (ocean-in) − (top-out) = Qᵇ − Qᵘ.
    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        T_top = @inbounds temperature[i, j, Nz]
        T_bottom = @inbounds temperature[i, j, 1]
        top = getflux(external_heat_fluxes.top, i, j, grid, T_top, clock, model_fields)
        bottom = getflux(external_heat_fluxes.bottom, i, j, grid, T_bottom, clock, model_fields)
        difference += bottom - top
    end

    return difference
end

function column_shortwave_flux_difference(thermodynamics::ColumnEnergyThermodynamics)
    compute_column_shortwave_flux!(thermodynamics)

    I = thermodynamics.auxiliary.shortwave_flux
    grid = thermodynamics.fields.internal_energy.grid
    Nz = size(grid, 3)
    difference = zero(eltype(grid))

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        difference += @inbounds I[i, j, Nz+1] - I[i, j, 1]
    end

    return difference
end

@inline function budget_relative_residual(residual, actual_change, expected_change)
    scale = max(abs(actual_change), abs(expected_change), one(abs(residual)))
    return abs(residual) / scale
end

"""
    column_energy_budget(thermodynamics, external_heat_fluxes, clock, model_fields, initial_energy, dt;
                         final_energy = column_integrated_energy(thermodynamics),
                         surface_stefan_residual_flux = 0)

Return a named tuple summarizing the column-integrated energy budget over one
step with duration `dt`. The boundary flux change is the requested (uncapped)
`external_heat_fluxes` difference; surface melt at a `MeltingConstrainedFluxBalance`
boundary is accounted separately through `surface_stefan_residual_flux`.
"""
function column_energy_budget(thermodynamics::ColumnEnergyThermodynamics,
                              external_heat_fluxes, clock, model_fields,
                              initial_energy,
                              Δt;
                              final_energy = column_integrated_energy(thermodynamics),
                              surface_stefan_residual_flux = 0)
    boundary_flux_change = Δt * column_boundary_energy_flux_difference(thermodynamics, external_heat_fluxes, clock, model_fields)
    shortwave_flux_change = Δt * column_shortwave_flux_difference(thermodynamics)
    surface_stefan_residual_change = Δt * surface_stefan_residual_flux
    expected_change = boundary_flux_change + shortwave_flux_change + surface_stefan_residual_change
    actual_change = final_energy - initial_energy
    residual = actual_change - expected_change

    return (initial = initial_energy,
            final = final_energy,
            actual_change,
            boundary_flux_change,
            shortwave_flux_change,
            surface_stefan_residual_change,
            expected_change,
            residual,
            relative_residual = budget_relative_residual(residual,
                                                         actual_change,
                                                         expected_change))
end

"""
    column_salt_budget(thermodynamics, initial_salt, dt; final_salt = column_integrated_salinity(thermodynamics))

Return a named tuple summarizing the closed-boundary column-integrated salt
budget over one step with duration `dt`.
"""
function column_salt_budget(thermodynamics::ColumnEnergyThermodynamics,
                            initial_salt,
                            Δt;
                            final_salt = column_integrated_salinity(thermodynamics))
    expected_change = zero(final_salt - initial_salt)
    actual_change = final_salt - initial_salt
    residual = actual_change - expected_change

    return (initial = initial_salt,
            final = final_salt,
            actual_change,
            expected_change,
            residual,
            relative_residual = budget_relative_residual(residual,
                                                         actual_change,
                                                         expected_change))
end
