@inline function store_tridiagonal_row!(auxiliary, i, j, k, Nz, left_factor, right_factor, boundary_factor)
    store_tridiagonal_row!(auxiliary, i, j, k, Nz, left_factor, left_factor, right_factor, right_factor, boundary_factor)
end

@inline function store_tridiagonal_row!(auxiliary, i, j, k, Nz, left_self, left_neighbor, right_self, right_neighbor, boundary_factor)
    @inbounds begin
        auxiliary.diagonal[i, j, k] = 1 + left_self + right_self + boundary_factor
        k > 1 && (auxiliary.lower_diagonal[i, j, k-1] = -left_neighbor)
        auxiliary.upper_diagonal[i, j, k] = -right_neighbor
        k == Nz && (auxiliary.lower_diagonal[i, j, k] = zero(left_self))
    end
end

@inline function column_energy_diffusion_factor(i, j, kc, kf, kn, grid, auxiliary, fields, κE, Δt)
    @inbounds D = auxiliary.thermal_conductivity[i, j, kf] * fields.temperature_energy_derivative[i, j, kn] + κE
    return diffusion_factor(i, j, kc, kf, grid, D, Δt)
end

@inline function column_affine_conductive_flux(i, j, kf, grid, auxiliary, fields)
    @inbounds begin
        kᵇ, kᵗ = kf - 1, kf
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
        if column_consolidated(consolidation_thickness, grid, i, j)
            κE = energy_diffusivity(energy_transport)
            I = auxiliary.shortwave_flux

            left_self     = k == 1  ? zero(Δt) : column_energy_diffusion_factor(i, j, k, k, k, grid, auxiliary, fields, κE, Δt)
            left_neighbor = k == 1  ? zero(Δt) : column_energy_diffusion_factor(i, j, k, k, k-1, grid, auxiliary, fields, κE, Δt)
            right_self     = k == Nz ? zero(Δt) : column_energy_diffusion_factor(i, j, k, k+1, k, grid, auxiliary, fields, κE, Δt)
            right_neighbor = k == Nz ? zero(Δt) : column_energy_diffusion_factor(i, j, k, k+1, k+1, grid, auxiliary, fields, κE, Δt)

            boundary_factor = zero(Δt)
            k == 1  && (boundary_factor += column_bottom_boundary_energy_factor(bcs.bottom, Qe.bottom, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt))
            k == Nz && (boundary_factor += column_top_boundary_energy_factor(bcs.top, Qe.top, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt))

            store_tridiagonal_row!(auxiliary, i, j, k, Nz, left_self, left_neighbor, right_self, right_neighbor, boundary_factor)

            left_flux = k == 1 ? column_bottom_boundary_energy_flux(bcs.bottom, Qe.bottom, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) :
                                 column_affine_conductive_flux(i, j, k, grid, auxiliary, fields)

            right_flux = k == Nz ? column_top_boundary_energy_flux(bcs.top, Qe.top, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) :
                                   column_affine_conductive_flux(i, j, k+1, grid, auxiliary, fields)

            flux_tendency = (right_flux + I[i, j, k+1] - left_flux - I[i, j, k]) / Δzᶜᶜᶜ(i, j, k, grid)

            # Conservative moving-grid balance: Jⁿ⁺¹ Δr Eⁿ⁺¹ = Jⁿ Δr Eⁿ + δzᵤ Eᵘᵖᵤ - δzₗ Eᵘᵖₗ + Δt (Fᵤ - Fₗ)
            moving_tendency = (moving_face_displacement_flux(i, j, k+1, grid, initial_energy, fields, relation, bcs.bottom) -
                               moving_face_displacement_flux(i, j, k, grid, initial_energy, fields, relation, bcs.bottom)) / Δzᶜᶜᶜ(i, j, k, grid)

            auxiliary.energy_rhs[i, j, k] = previous_to_current_column_metric_ratio(i, j, k, grid) * initial_energy[i, j, k] + moving_tendency + Δt * flux_tendency
        else
            # Unconsolidated columns hold their internal energy; the slab balance in the volume update takes the growth.
            store_tridiagonal_row!(auxiliary, i, j, k, Nz, zero(Δt), zero(Δt), zero(Δt))
            auxiliary.energy_rhs[i, j, k] = initial_energy[i, j, k]
        end
    end
end

"""
    assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)

Assemble the tridiagonal backward-Euler system of one internal-energy step. Columns thinner than
`consolidation_thickness` get an identity row.
"""
function assemble_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt; initial_energy = thermodynamics.fields.internal_energy)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xyz, _assemble_column_energy_system!,
            thermodynamics.auxiliary, thermodynamics.fields, initial_energy, grid, thermodynamics.heat_boundary_conditions,
            external_heat_fluxes, thermodynamics.energy_transport, thermodynamics.relation, consolidation_thickness, clock, model_fields, Δt)
    return nothing
end

#####
##### Surface temperature of a `MeltingConstrainedFluxBalance` top: the massless-surface balance
##### Qᵘ(Tₛ) = conductance (T_top − Tₛ), capped at melting.
#####

# A Float32 flux difference cannot be resolved to an absolute 1e-6 W/m²
# when either opposing flux is hundreds of W/m². Scale the root tolerance
# with the roundoff of the two fluxes being balanced.
@inline function column_flux_balance_tolerance(grid, external_flux, conductive_flux, conductance, temperature_scale)
    FT = eltype(grid)
    scale = max(one(external_flux), abs(external_flux), abs(conductive_flux))
    return max(oftype(scale, 1e-6), 16eps(FT) * scale,
               8eps(FT) * abs(conductance) * max(one(temperature_scale), abs(temperature_scale)))
end

@inline function column_surface_temperature_balance(i, j, grid, Tₛ⁻, conductance, T_top, external_flux, clock, model_fields)
    FT = eltype(grid)
    conductance > eps(FT) || return Tₛ⁻
    balance(T) = getflux(external_flux, i, j, grid, T, clock, model_fields) - conductance * (T_top - T)
    tolerance(T) = column_flux_balance_tolerance(grid,
                                                  getflux(external_flux, i, j, grid, T, clock, model_fields),
                                                  conductance * (T_top - T), conductance,
                                                  max(abs(T_top), abs(T)))

    # Keep the fast solve when it actually closes the heat balance. For thin
    # snow, secant can report nonconvergence and return a remote iterate even
    # when the physical root lies close to melting.
    solution = find_zero(balance, SecantMethod{FT}(Tₛ⁻ + one(FT), Tₛ⁻), CompactSolution())
    candidate = solution.root
    if isfinite(candidate)
        residual = balance(candidate)
        isfinite(residual) && abs(residual) <= tolerance(candidate) && return candidate
    end

    # The surface is capped at the melting point by the caller. A root above
    # 0°C therefore needs no solve, while a colder root can be bracketed.
    warm = zero(FT)
    warm_residual = balance(warm)
    isfinite(warm_residual) || return FT(NaN)
    warm_residual <= 0 && return warm

    cold = min(T_top, warm) - one(FT)
    width = one(FT)
    cold_residual = balance(cold)
    for _ in 1:32
        isfinite(cold_residual) || return FT(NaN)
        cold_residual <= 0 && break
        width *= 2
        cold -= width
        cold_residual = balance(cold)
    end
    (isfinite(cold_residual) && cold_residual <= 0) || return FT(NaN)

    for _ in 1:64
        trial = (cold + warm) / 2
        residual = balance(trial)
        isfinite(residual) || return FT(NaN)
        abs(residual) <= tolerance(trial) && return trial
        if residual > 0
            warm = trial
        else
            cold = trial
        end
    end

    trial = (cold + warm) / 2
    return abs(balance(trial)) <= tolerance(trial) ? trial : FT(NaN)
end

@kernel function _update_column_surface_temperature!(auxiliary, fields, grid, external_flux, relation, clock, model_fields)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        # Masked cells and zero-height columns have no ice surface balance.
        # Keep a finite diagnostic value without solving an undefined balance.
        if Oceananigans.ImmersedBoundaries.immersed_cell(i, j, Nz, grid) ||
           column_height(grid, i, j) <= eps(eltype(grid))
            auxiliary.surface_temperature[i, j, 1] = fields.temperature[i, j, Nz]
        else
            conductance = column_boundary_temperature_conductance(i, j, Nz+1, Nz, grid, auxiliary)
            Tₘ = melting_temperature(relation.phase_transitions.liquidus, fields.bulk_salinity[i, j, Nz])
            Tₛ = column_surface_temperature_balance(i, j, grid, auxiliary.surface_temperature[i, j, 1], conductance,
                                                    fields.temperature[i, j, Nz], external_flux, clock, model_fields)
            auxiliary.surface_temperature[i, j, 1] = min(Tₛ, Tₘ)
        end
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
