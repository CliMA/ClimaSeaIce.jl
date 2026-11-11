# Conduction through the column's bottom and top faces comes from the boundary conditions instead, so the interior
# terms vanish there, as in Oceananigans' implicit vertical diffusion.
@inline interior_column_face(kf, grid, value) = ifelse((kf > 1) & (kf <= size(grid, 3)), value, zero(value))

# Row k of a tridiagonal system, or an identity row that keeps `fixed_value` in an inactive cell. `lower` is row k+1's
# coefficient of cell k, as in Oceananigans' BatchedTridiagonalSolver.
@inline function store_tridiagonal_row!(i, j, k, auxiliary, active, diagonal, lower, upper, rhs, fixed_value)
    d, l, u, r = ifelse(active, (diagonal, -lower, -upper, rhs), (one(diagonal), zero(lower), zero(upper), fixed_value))
    @inbounds begin
        auxiliary.diagonal[i, j, k] = d
        auxiliary.lower_diagonal[i, j, k] = l
        auxiliary.upper_diagonal[i, j, k] = u
        auxiliary.energy_rhs[i, j, k] = r
    end
    return nothing
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

# Heat flux through face kf: the boundary conditions on the column's bottom and top faces, conduction in between.
@inline function column_face_energy_flux(i, j, kf, grid, bcs, Qe, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    bottom = column_bottom_boundary_energy_flux(i, j, 1, grid, bcs.bottom, Qe.bottom, auxiliary, fields, relation, clock, model_fields, Δt)
    top = column_top_boundary_energy_flux(i, j, Nz, grid, bcs.top, Qe.top, auxiliary, fields, relation, clock, model_fields, Δt)
    return ifelse(kf == 1, bottom, ifelse(kf == Nz + 1, top, column_affine_conductive_flux(i, j, kf, grid, auxiliary, fields)))
end

@kernel function _assemble_column_energy_system!(auxiliary, fields, initial_energy, grid, heat_boundary_conditions, external_heat_fluxes, energy_transport, relation, consolidation_thickness, clock, model_fields, Δt)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)
    bcs = heat_boundary_conditions
    Qe = external_heat_fluxes

    @inbounds begin
        κE = energy_diffusivity(energy_transport)
        I = auxiliary.shortwave_flux
        kᵘ = min(Nz, k+1)

        # Implicit conduction through the bottom face k and the top face k+1 of cell k: the factors of its own energy,
        # and of the neighbours' energies in the off-diagonals.
        bottom_self = interior_column_face(k,   grid, column_energy_diffusion_factor(i, j, k,  k,   k,  grid, auxiliary, fields, κE, Δt))
        top_self    = interior_column_face(k+1, grid, column_energy_diffusion_factor(i, j, k,  k+1, k,  grid, auxiliary, fields, κE, Δt))
        upper       = interior_column_face(k+1, grid, column_energy_diffusion_factor(i, j, k,  k+1, kᵘ, grid, auxiliary, fields, κE, Δt))
        lower       = interior_column_face(k+1, grid, column_energy_diffusion_factor(i, j, kᵘ, k+1, k,  grid, auxiliary, fields, κE, Δt))
        boundary_factor = ifelse(k == 1, column_bottom_boundary_energy_factor(i, j, k, grid, bcs.bottom, Qe.bottom, auxiliary, fields, relation, clock, model_fields, Δt), zero(Δt)) +
                          ifelse(k == Nz, column_top_boundary_energy_factor(i, j, k, grid, bcs.top, Qe.top, auxiliary, fields, relation, clock, model_fields, Δt), zero(Δt))
        diagonal = 1 + bottom_self + top_self + boundary_factor

        Fᵇ = column_face_energy_flux(i, j, k,   grid, bcs, Qe, auxiliary, fields, relation, clock, model_fields, Δt)
        Fᵗ = column_face_energy_flux(i, j, k+1, grid, bcs, Qe, auxiliary, fields, relation, clock, model_fields, Δt)
        flux_tendency = (Fᵗ + I[i, j, k+1] - Fᵇ - I[i, j, k]) / Δzᶜᶜᶜ(i, j, k, grid)

        # Conservative moving-grid balance: Jⁿ⁺¹ Δr Eⁿ⁺¹ = Jⁿ Δr Eⁿ + δzᵤ Eᵘᵖᵤ - δzₗ Eᵘᵖₗ + Δt (Fᵤ - Fₗ)
        moving_tendency = (moving_face_displacement_flux(i, j, k+1, grid, initial_energy, fields, relation, bcs.bottom, model_fields) -
                           moving_face_displacement_flux(i, j, k, grid, initial_energy, fields, relation, bcs.bottom, model_fields)) / Δzᶜᶜᶜ(i, j, k, grid)
        rhs = previous_to_current_column_metric_ratio(i, j, k, grid) * initial_energy[i, j, k] + moving_tendency + Δt * flux_tendency

        # Columns thinner than the consolidation thickness take the basal-reservoir state, capped at melting.
        salinity = fields.bulk_salinity[i, j, k]
        melting_point = melting_temperature(relation.phase_transitions.liquidus, salinity)
        reservoir_temperature = min(column_bottom_reference_temperature(i, j, grid, bcs.bottom, fields, relation, model_fields), melting_point)
        reservoir_energy = internal_energy(relation, reservoir_temperature, salinity)
        consolidated = column_consolidated(column_height(i, j, grid), consolidation_thickness[i, j, 1])
        store_tridiagonal_row!(i, j, k, auxiliary, consolidated, diagonal, lower, upper, rhs, reservoir_energy)
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
    kᵘ = min(Nz, k+1)

    bottom = interior_column_face(k,   grid, diffusion_factor(i, j, k,  k,   grid, κS, Δt))
    top    = interior_column_face(k+1, grid, diffusion_factor(i, j, k,  k+1, grid, κS, Δt))
    lower  = interior_column_face(k+1, grid, diffusion_factor(i, j, kᵘ, k+1, grid, κS, Δt))

    moving_tendency = (moving_face_displacement_flux(i, j, k+1, grid, S, fields, nothing, nothing, nothing) -
                       moving_face_displacement_flux(i, j, k, grid, S, fields, nothing, nothing, nothing)) / Δzᶜᶜᶜ(i, j, k, grid)
    @inbounds rhs = previous_to_current_column_metric_ratio(i, j, k, grid) * S[i, j, k] + moving_tendency
    active = column_height(i, j, grid) > 0
    @inbounds store_tridiagonal_row!(i, j, k, auxiliary, active, 1 + bottom + top, lower, top, rhs, S[i, j, k])
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
