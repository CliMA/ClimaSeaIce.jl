@inline function store_tridiagonal_row!(auxiliary, i, j, k, Nz, left_factor, right_factor, boundary_factor)
    @inbounds begin
        auxiliary.diagonal[i, j, k] = 1 + left_factor + right_factor + boundary_factor
        k > 1 && (auxiliary.lower_diagonal[i, j, k-1] = -left_factor)
        auxiliary.upper_diagonal[i, j, k] = -right_factor
        k == Nz && (auxiliary.lower_diagonal[i, j, k] = zero(left_factor))
    end
end

@kernel function _assemble_column_energy_system!(auxiliary, fields, grid, heat_boundary_conditions, external_heat_fluxes, relation, consolidation_thickness, clock, model_fields, Δt)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)
    bcs = heat_boundary_conditions
    Qe = external_heat_fluxes

    @inbounds begin
        if column_consolidated(consolidation_thickness, grid, i, j)
            D = auxiliary.effective_energy_diffusivity
            C = auxiliary.salinity_coupling_diffusivity
            I = auxiliary.shortwave_flux
            S = fields.bulk_salinity
            E = fields.internal_energy

            left_factor  = k == 1  ? zero(Δt) : diffusion_factor(i, j, k, k, grid, D, Δt)
            right_factor = k == Nz ? zero(Δt) : diffusion_factor(i, j, k, k+1, grid, D, Δt)

            boundary_factor = zero(Δt)
            k == 1  && (boundary_factor += column_bottom_boundary_energy_factor(bcs.bottom, Qe.bottom, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt))
            k == Nz && (boundary_factor += column_top_boundary_energy_factor(bcs.top, Qe.top, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt))

            store_tridiagonal_row!(auxiliary, i, j, k, Nz, left_factor, right_factor, boundary_factor)

            left_flux = k == 1 ? column_bottom_boundary_energy_flux(bcs.bottom, Qe.bottom, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) :
                                 salinity_coupling_flux(i, j, k, grid, C, S)

            right_flux = k == Nz ? column_top_boundary_energy_flux(bcs.top, Qe.top, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) :
                                   salinity_coupling_flux(i, j, k+1, grid, C, S)

            flux_tendency = (right_flux + I[i, j, k+1] - left_flux - I[i, j, k]) / Δzᶜᶜᶜ(i, j, k, grid)

            # Conservative moving-grid balance: Jⁿ⁺¹ Δr Eⁿ⁺¹ = Jⁿ Δr Eⁿ + δzᵤ Eᵘᵖᵤ - δzₗ Eᵘᵖₗ + Δt (Fᵤ - Fₗ)
            moving_tendency = (moving_face_displacement_flux(i, j, k+1, grid, E, fields, relation, bcs.bottom) -
                               moving_face_displacement_flux(i, j, k, grid, E, fields, relation, bcs.bottom)) / Δzᶜᶜᶜ(i, j, k, grid)

            auxiliary.energy_rhs[i, j, k] = previous_to_current_column_metric_ratio(i, j, k, grid) * E[i, j, k] + moving_tendency + Δt * flux_tendency
        else
            # Unconsolidated columns hold their internal energy; the slab balance in the volume update takes the growth.
            store_tridiagonal_row!(auxiliary, i, j, k, Nz, zero(Δt), zero(Δt), zero(Δt))
            auxiliary.energy_rhs[i, j, k] = fields.internal_energy[i, j, k]
        end
    end
end

"""
    assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)

Assemble the tridiagonal backward-Euler system of one internal-energy step. Columns thinner than
`consolidation_thickness` get an identity row.
"""
function assemble_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xyz, _assemble_column_energy_system!,
            thermodynamics.auxiliary, thermodynamics.fields, grid, thermodynamics.heat_boundary_conditions,
            external_heat_fluxes, thermodynamics.relation, consolidation_thickness, clock, model_fields, Δt)
    return nothing
end

#####
##### Surface temperature of a `MeltingConstrainedFluxBalance` top: the massless-surface balance
##### Qᵘ(Tₛ) = conductance (T_top − Tₛ), capped at melting.
#####

@inline function column_surface_temperature_balance(i, j, grid, Tₛ⁻, conductance, T_top, external_flux, clock, model_fields)
    FT = eltype(grid)
    conductance > eps(FT) || return Tₛ⁻
    balance(T) = getflux(external_flux, i, j, grid, T, clock, model_fields) - conductance * (T_top - T)
    return find_zero(balance, SecantMethod{FT}(Tₛ⁻ + one(FT), Tₛ⁻), CompactSolution()).root
end

@kernel function _update_column_surface_temperature!(auxiliary, fields, grid, external_flux, relation, clock, model_fields)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        conductance = column_boundary_temperature_conductance(i, j, Nz+1, Nz, grid, auxiliary)
        Tₘ = melting_temperature(relation.phase_transitions.liquidus, fields.bulk_salinity[i, j, Nz])
        Tₛ = column_surface_temperature_balance(i, j, grid, auxiliary.surface_temperature[i, j, 1], conductance,
                                                fields.temperature[i, j, Nz], external_flux, clock, model_fields)
        auxiliary.surface_temperature[i, j, 1] = min(Tₛ, Tₘ)
    end
end

"""
    solve_column_energy_system!(thermodynamics)

Solve the assembled column energy system into the internal-energy field.
"""
solve_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics) = solve!(thermodynamics.fields.internal_energy, thermodynamics.solvers.energy_solver, thermodynamics.auxiliary.energy_rhs)

@kernel function _assemble_column_salinity_system!(auxiliary, fields, grid, κS, Δt)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)
    S = fields.bulk_salinity

    left_factor  = k == 1  ? zero(Δt) : diffusion_factor(i, j, k, k, grid, κS, Δt)
    right_factor = k == Nz ? zero(Δt) : diffusion_factor(i, j, k, k+1, grid, κS, Δt)
    store_tridiagonal_row!(auxiliary, i, j, k, Nz, left_factor, right_factor, zero(Δt))

    moving_tendency = (moving_face_displacement_flux(i, j, k+1, grid, S, fields, nothing, nothing) -
                       moving_face_displacement_flux(i, j, k, grid, S, fields, nothing, nothing)) / Δzᶜᶜᶜ(i, j, k, grid)

    @inbounds auxiliary.energy_rhs[i, j, k] = previous_to_current_column_metric_ratio(i, j, k, grid) * S[i, j, k] + moving_tendency
end

"""
    assemble_column_salinity_system!(thermodynamics, Δt)

Assemble the tridiagonal backward-Euler system of one closed-boundary bulk-salinity diffusion step.
"""
function assemble_column_salinity_system!(thermodynamics::ColumnEnergyThermodynamics, Δt)
    grid = thermodynamics.fields.bulk_salinity.grid
    κS = salinity_diffusivity(thermodynamics.salinity_transport)
    launch!(architecture(grid), grid, :xyz, _assemble_column_salinity_system!, thermodynamics.auxiliary, thermodynamics.fields, grid, κS, Δt)
    return nothing
end

"""
    solve_column_salinity_system!(thermodynamics)

Solve the assembled bulk-salinity system into the bulk-salinity field.
"""
solve_column_salinity_system!(thermodynamics::ColumnEnergyThermodynamics) = solve!(thermodynamics.fields.bulk_salinity, thermodynamics.solvers.salinity_solver, thermodynamics.auxiliary.energy_rhs)

"""
    column_salinity_time_step!(thermodynamics, Δt)

Advance the prognostic bulk salinity by one conservative moving-grid diffusion step. Prescribed salinity is left unchanged.
"""
column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics, Δt) = column_salinity_time_step!(thermodynamics, thermodynamics.salinity_closure, Δt)
column_salinity_time_step!(thermodynamics, ::PrescribedBulkSalinity, Δt) = nothing

function column_salinity_time_step!(thermodynamics, ::PrognosticBulkSalinity, Δt)
    assemble_column_salinity_system!(thermodynamics, Δt)
    solve_column_salinity_system!(thermodynamics)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    return nothing
end
