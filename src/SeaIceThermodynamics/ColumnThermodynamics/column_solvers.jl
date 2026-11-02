@inline function store_tridiagonal_row!(i, j, k, grid, auxiliary, active, left_self, right_self, lower_neighbor, upper_neighbor, boundary_factor)
    @inbounds begin
        diagonal = 1 + left_self + right_self + boundary_factor
        auxiliary.diagonal[i, j, k] = ifelse(active, diagonal, one(diagonal))
        auxiliary.lower_diagonal[i, j, k] = ifelse(active, -lower_neighbor, zero(lower_neighbor))
        auxiliary.upper_diagonal[i, j, k] = ifelse(active, -upper_neighbor, zero(upper_neighbor))
    end
end

@inline function column_energy_diffusion_factor(i, j, kc, kf, kn, grid, auxiliary, fields, κE, Δt)
    @inbounds D = auxiliary.thermal_conductivity[i, j, kf] * fields.temperature_energy_derivative[i, j, kn] + κE
    return diffusion_factor(i, j, kc, kf, grid, D, Δt)
end

@inline function column_affine_conductive_flux(i, j, kf, grid, auxiliary, fields)
    @inbounds begin
        kᵇ, kᵗ = max(kf - 1, 1), min(kf, size(grid, 3))
        T = fields.temperature
        E = fields.internal_energy
        TE = fields.temperature_energy_derivative
        K = auxiliary.thermal_conductivity[i, j, kf]
        return K * ((T[i, j, kᵗ] - TE[i, j, kᵗ] * E[i, j, kᵗ]) -
                    (T[i, j, kᵇ] - TE[i, j, kᵇ] * E[i, j, kᵇ])) / Δzᶜᶜᶠ(i, j, kf, grid)
    end
end

@kernel function _assemble_column_energy_system!(auxiliary, fields, initial_energy, grid, heat_boundary_conditions, external_heat_fluxes, energy_transport, relation, consolidation_thickness, clock, model_fields, Δt)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)
    bcs = heat_boundary_conditions
    Qe = external_heat_fluxes

    @inbounds begin
        κE = energy_diffusivity(energy_transport)
        I = auxiliary.shortwave_flux

        left_self     = ifelse(k == 1, zero(Δt), column_energy_diffusion_factor(i, j, k, k, k, grid, auxiliary, fields, κE, Δt))
        lower_neighbor = ifelse(k == Nz, zero(Δt), column_energy_diffusion_factor(i, j, min(Nz, k+1), k+1, k, grid, auxiliary, fields, κE, Δt))
        right_self     = ifelse(k == Nz, zero(Δt), column_energy_diffusion_factor(i, j, k, k+1, k, grid, auxiliary, fields, κE, Δt))
        right_neighbor = ifelse(k == Nz, zero(Δt), column_energy_diffusion_factor(i, j, k, k+1, min(Nz, k+1), grid, auxiliary, fields, κE, Δt))

        boundary_factor = ifelse(k == 1, column_bottom_boundary_energy_factor(i, j, k, grid, bcs.bottom, Qe.bottom, auxiliary, fields, relation, clock, model_fields, Δt), zero(Δt)) +
                          ifelse(k == Nz, column_top_boundary_energy_factor(i, j, k, grid, bcs.top, Qe.top, auxiliary, fields, relation, clock, model_fields, Δt), zero(Δt))
        consolidated = column_consolidated(column_height(grid, i, j), consolidation_thickness[i, j, 1])

        store_tridiagonal_row!(i, j, k, grid, auxiliary, consolidated, left_self, right_self, lower_neighbor, right_neighbor, boundary_factor)

        left_flux = ifelse(k == 1, column_bottom_boundary_energy_flux(i, j, k, grid, bcs.bottom, Qe.bottom, auxiliary, fields, relation, clock, model_fields, Δt),
                          column_affine_conductive_flux(i, j, k, grid, auxiliary, fields))

        right_flux = ifelse(k == Nz, column_top_boundary_energy_flux(i, j, k, grid, bcs.top, Qe.top, auxiliary, fields, relation, clock, model_fields, Δt),
                           column_affine_conductive_flux(i, j, k+1, grid, auxiliary, fields))

        flux_tendency = (right_flux + I[i, j, k+1] - left_flux - I[i, j, k]) / Δzᶜᶜᶜ(i, j, k, grid)

        # Conservative moving-grid balance: Jⁿ⁺¹ Δr Eⁿ⁺¹ = Jⁿ Δr Eⁿ + δzᵤ Eᵘᵖᵤ - δzₗ Eᵘᵖₗ + Δt (Fᵤ - Fₗ)
        moving_tendency = (moving_face_displacement_flux(i, j, k+1, grid, initial_energy, fields, relation, bcs.bottom) -
                           moving_face_displacement_flux(i, j, k, grid, initial_energy, fields, relation, bcs.bottom)) / Δzᶜᶜᶜ(i, j, k, grid)

        salinity = fields.bulk_salinity[i, j, k]
        melting_point = melting_temperature(relation.phase_transitions.liquidus, salinity)
        reservoir_temperature = min(column_bottom_reference_temperature(bcs.bottom, i, j, grid, fields, relation), melting_point)
        reservoir_energy = internal_energy(relation, reservoir_temperature, salinity)
        auxiliary.energy_rhs[i, j, k] = ifelse(consolidated, previous_to_current_column_metric_ratio(i, j, k, grid) * initial_energy[i, j, k] + moving_tendency + Δt * flux_tendency, reservoir_energy)
    end
end

"""
    assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)

Assemble the tridiagonal backward-Euler system of one internal-energy step. Columns thinner than
`consolidation_thickness` take the basal-reservoir state capped at the local melting temperature.
"""
function assemble_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt; initial_energy = thermodynamics.fields.internal_energy)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xyz, _assemble_column_energy_system!,
            thermodynamics.auxiliary, thermodynamics.fields, initial_energy, grid, thermodynamics.heat_boundary_conditions,
            external_heat_fluxes, thermodynamics.energy_transport, thermodynamics.relation, consolidation_thickness, clock, model_fields, Δt)
    return nothing
end

@kernel function _column_energy_residual!(auxiliary, E, initial_energy, grid, relative_tolerance)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)
    @inbounds begin
        lower = ifelse(k > 1, auxiliary.lower_diagonal[i, j, max(1, k-1)] * E[i, j, max(1, k-1)], zero(eltype(grid)))
        upper = ifelse(k < Nz, auxiliary.upper_diagonal[i, j, k] * E[i, j, min(Nz, k+1)], zero(eltype(grid)))
        diagonal = auxiliary.diagonal[i, j, k] * E[i, j, k]
        rhs = auxiliary.energy_rhs[i, j, k]
        scale = max(1, abs(E[i, j, k] - initial_energy[i, j, k]))
        tolerance = max(relative_tolerance * scale,
                        4eps(eltype(grid)) * (abs(diagonal) + abs(lower) + abs(upper) + abs(rhs)))
        residual = abs(diagonal + lower + upper - rhs) / tolerance
        auxiliary.energy_rhs[i, j, k] = ifelse(Oceananigans.ImmersedBoundaries.immersed_cell(i, j, k, grid), zero(residual), residual)
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

    left_factor  = ifelse(k == 1, zero(Δt), diffusion_factor(i, j, k, k, grid, κS, Δt))
    right_factor = ifelse(k == Nz, zero(Δt), diffusion_factor(i, j, k, k+1, grid, κS, Δt))
    lower_factor = ifelse(k == Nz, zero(Δt), diffusion_factor(i, j, min(Nz, k+1), k+1, grid, κS, Δt))
    active = column_height(grid, i, j) > 0
    store_tridiagonal_row!(i, j, k, grid, auxiliary, active, left_factor, right_factor, lower_factor, right_factor, zero(Δt))

    moving_tendency = (moving_face_displacement_flux(i, j, k+1, grid, S, fields, nothing, nothing) -
                       moving_face_displacement_flux(i, j, k, grid, S, fields, nothing, nothing)) / Δzᶜᶜᶜ(i, j, k, grid)

    @inbounds auxiliary.energy_rhs[i, j, k] = ifelse(active, previous_to_current_column_metric_ratio(i, j, k, grid) * S[i, j, k] + moving_tendency, S[i, j, k])
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
function column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics, Δt)
    thermodynamics.salinity_closure isa PrognosticBulkSalinity || return nothing
    assemble_column_salinity_system!(thermodynamics, Δt)
    solve_column_salinity_system!(thermodynamics)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    return nothing
end
