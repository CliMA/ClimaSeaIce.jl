@kernel function _assemble_column_energy_system!(auxiliary, fields, grid, heat_boundary_conditions, external_heat_fluxes, relation, consolidation_thickness, clock, model_fields, Δt)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)

    if column_consolidated(consolidation_thickness, grid, i, j)
    @inbounds begin
        D = auxiliary.effective_energy_diffusivity
        C = auxiliary.salinity_coupling_diffusivity
        I = auxiliary.shortwave_flux
        S = fields.bulk_salinity

        left_factor = if k == 1
            zero(Δt)
        else
            diffusion_factor(i, j, k, k, grid, D, Δt)
        end

        right_factor = if k == Nz
            zero(Δt)
        else
            diffusion_factor(i, j, k, k+1, grid, D, Δt)
        end

        bottom_factor = if k == 1
            column_bottom_boundary_energy_factor(heat_boundary_conditions.bottom, external_heat_fluxes.bottom,
                                                 i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
        else
            zero(Δt)
        end

        top_factor = if k == Nz
            column_top_boundary_energy_factor(heat_boundary_conditions.top, external_heat_fluxes.top,
                                              i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
        else
            zero(Δt)
        end

        auxiliary.diagonal[i, j, k] = 1 + left_factor + right_factor + bottom_factor + top_factor

        if k > 1
            auxiliary.lower_diagonal[i, j, k-1] = -left_factor
        end

        if k < Nz
            auxiliary.upper_diagonal[i, j, k] = -right_factor
        else
            auxiliary.upper_diagonal[i, j, k] = zero(right_factor)
            auxiliary.lower_diagonal[i, j, k] = zero(left_factor)
        end

        left_flux = if k == 1
            column_bottom_boundary_energy_flux(heat_boundary_conditions.bottom, external_heat_fluxes.bottom,
                                               i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
        else
            salinity_coupling_flux(i, j, k, grid, C, S)
        end
        left_flux += I[i, j, k]

        right_flux = if k == Nz
            column_top_boundary_energy_flux(heat_boundary_conditions.top, external_heat_fluxes.top,
                                            i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
        else
            salinity_coupling_flux(i, j, k+1, grid, C, S)
        end
        right_flux += I[i, j, k+1]

        flux_tendency = (right_flux - left_flux) / Δzᶜᶜᶜ(i, j, k, grid)

        # Conservative moving-grid balance divided by the current layer thickness:
        # Jⁿ⁺¹ Δr Eⁿ⁺¹ = Jⁿ Δr Eⁿ + δzᵤ Eᵘᵖᵤ - δzₗ Eᵘᵖₗ + Δt (Fᵤ - Fₗ).
        moving_flux_difference =
            moving_face_displacement_flux(i, j, k+1, grid, fields.internal_energy, fields, relation,
                                          heat_boundary_conditions.bottom,
                                          heat_boundary_conditions.top) -
            moving_face_displacement_flux(i, j, k, grid, fields.internal_energy, fields, relation,
                                          heat_boundary_conditions.bottom,
                                          heat_boundary_conditions.top)
        moving_tendency = moving_flux_difference / Δzᶜᶜᶜ(i, j, k, grid)
        metric_ratio = previous_to_current_column_metric_ratio(i, j, k, grid)
        auxiliary.energy_rhs[i, j, k] =
            metric_ratio * fields.internal_energy[i, j, k] + moving_tendency + Δt * flux_tendency
    end
    else
    @inbounds begin
        # Unconsolidated column: no resolved interior profile. Hold the internal energy with an identity row;
        # the thickness change is taken up by the slab flux balance in the volume update.
        auxiliary.diagonal[i, j, k] = one(Δt)
        if k > 1
            auxiliary.lower_diagonal[i, j, k-1] = zero(Δt)
        end
        if k < Nz
            auxiliary.upper_diagonal[i, j, k] = zero(Δt)
        else
            auxiliary.lower_diagonal[i, j, k] = zero(Δt)
        end
        auxiliary.energy_rhs[i, j, k] = fields.internal_energy[i, j, k]
    end
    end
end

"""
    assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, dt)

Assemble the tridiagonal backward-Euler system for one internal-energy step.
`external_heat_fluxes` is the model's `(top, bottom)` flux set, evaluated through
`getflux` with `clock` and `model_fields`. A `consolidation_thickness` of `nothing`
treats every column as consolidated (the standalone path); otherwise columns thinner
than it get an identity row (no resolved profile).
"""
function assemble_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics,
                                        external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _assemble_column_energy_system!,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            grid,
            thermodynamics.heat_boundary_conditions,
            external_heat_fluxes,
            thermodynamics.relation,
            consolidation_thickness,
            clock,
            model_fields,
            Δt)

    return nothing
end

#####
##### Implicit surface-temperature solve for a `MeltingConstrainedFluxBalance` top, Icepack `temperature_changes`
##### style: an outer iteration around the tridiagonal solve makes the lagged top flux self-consistent with the
##### surface temperature `Tₛ`, found from the massless-surface balance `Q_atm(Tₛ) + conductance·(T_top − Tₛ) = 0`
##### and capped at melting. Negligible conduction decouples the surface, so `Tₛ` is left unchanged there.
#####

@inline function column_surface_temperature_balance(i, j, grid, Tₛ⁻, conductance, T_top, external_flux, clock, model_fields)
    FT = eltype(grid)
    conductance > eps(FT) || return Tₛ⁻
    # Massless-surface balance in the upward-positive convention: Qᵘ(Tₛ) = conductance·(T_top − Tₛ).
    balance(T) = getflux(external_flux, i, j, grid, T, clock, model_fields) - conductance * (T_top - T)
    solution = find_zero(balance, SecantMethod{FT}(Tₛ⁻ + one(FT), Tₛ⁻), CompactSolution())
    return solution.root
end

@kernel function _update_column_surface_temperature!(auxiliary, fields, grid, external_flux, relation, clock, model_fields)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        Tₛ⁻ = auxiliary.surface_temperature[i, j, 1]
        T_top = fields.temperature[i, j, Nz]
        S_top = fields.bulk_salinity[i, j, Nz]
        conductance = column_boundary_temperature_conductance(i, j, Nz+1, Nz, grid, auxiliary)
        Tₘ = melting_temperature(relation.phase_transitions.liquidus, S_top)
        Tₛ = column_surface_temperature_balance(i, j, grid, Tₛ⁻, conductance, T_top, external_flux, clock, model_fields)
        auxiliary.surface_temperature[i, j, 1] = min(Tₛ, Tₘ)
    end
end

"""
    solve_column_energy_system!(thermodynamics)

Solve the assembled column energy system into the internal-energy field.
"""
function solve_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics)
    solve!(thermodynamics.fields.internal_energy,
           thermodynamics.solvers.energy_solver,
           thermodynamics.auxiliary.energy_rhs)

    return nothing
end

@kernel function _assemble_column_salinity_system!(auxiliary, fields, grid, Δt)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        K = auxiliary.salinity_diffusivity

        left_factor = if k == 1
            zero(Δt)
        else
            diffusion_factor(i, j, k, k, grid, K, Δt)
        end

        right_factor = if k == Nz
            zero(Δt)
        else
            diffusion_factor(i, j, k, k+1, grid, K, Δt)
        end

        auxiliary.diagonal[i, j, k] = 1 + left_factor + right_factor

        if k > 1
            auxiliary.lower_diagonal[i, j, k-1] = -left_factor
        end

        if k < Nz
            auxiliary.upper_diagonal[i, j, k] = -right_factor
        else
            auxiliary.upper_diagonal[i, j, k] = zero(right_factor)
            auxiliary.lower_diagonal[i, j, k] = zero(left_factor)
        end

        metric_ratio = previous_to_current_column_metric_ratio(i, j, k, grid)

        # Same conservative moving-grid update as energy, without boundary heat
        # fluxes: Jⁿ⁺¹ Δr Sⁿ⁺¹ = Jⁿ Δr Sⁿ + δzᵤ Sᵘᵖᵤ - δzₗ Sᵘᵖₗ.
        moving_flux_difference =
            moving_face_displacement_flux(i, j, k+1, grid, fields.bulk_salinity) -
            moving_face_displacement_flux(i, j, k, grid, fields.bulk_salinity)
        moving_tendency = moving_flux_difference / Δzᶜᶜᶜ(i, j, k, grid)
        auxiliary.energy_rhs[i, j, k] = metric_ratio * fields.bulk_salinity[i, j, k] + moving_tendency
    end
end

"""
    assemble_column_salinity_system!(thermodynamics, dt)

Assemble the tridiagonal backward-Euler system for one closed-boundary bulk
salinity diffusion step.
"""
function assemble_column_salinity_system!(thermodynamics::ColumnEnergyThermodynamics, Δt)
    grid = thermodynamics.fields.bulk_salinity.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _assemble_column_salinity_system!,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            grid,
            Δt)

    return nothing
end

"""
    solve_column_salinity_system!(thermodynamics)

Solve the assembled scalar salinity system into the bulk-salinity field.
"""
function solve_column_salinity_system!(thermodynamics::ColumnEnergyThermodynamics)
    solve!(thermodynamics.fields.bulk_salinity,
           thermodynamics.solvers.salinity_solver,
           thermodynamics.auxiliary.energy_rhs)

    return nothing
end

"""
    column_salinity_time_step!(thermodynamics, dt)

Advance prognostic bulk salinity by one scalar diffusion step. On a
`MutableVerticalDiscretization`, the step evolves layer-integrated salinity
with the moving vertical metric. The step is a no-op for prescribed salinity.
"""
column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics, Δt) =
    column_salinity_time_step!(thermodynamics,
                               thermodynamics.salinity_closure,
                               thermodynamics.salinity_transport,
                               Δt)

function column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics,
                                    ::PrescribedBulkSalinity,
                                    salinity_transport,
                                    Δt)
    return nothing
end

function column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics,
                                    ::PrognosticBulkSalinity,
                                    ::BrineSalinityDiffusion,
                                    Δt)
    throw(ArgumentError("BrineSalinityDiffusion is a marker for future brine-salinity transport; use BulkSalinityDiffusion for the implemented scalar bulk-salinity step."))
end

# Both transports take the same closed-boundary moving-grid solve: `NoSalinityTransport` has zero diffusivity, so
# the tridiagonal system reduces to the identity and the step is the conservative metric remap alone, while
# `BulkSalinityDiffusion` adds the diffusion factors. No separate metric-remap path is needed.
function column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics,
                                    ::PrognosticBulkSalinity,
                                    salinity_transport,
                                    Δt)
    compute_column_salinity_diffusivity!(thermodynamics)
    assemble_column_salinity_system!(thermodynamics, Δt)
    solve_column_salinity_system!(thermodynamics)
    compute_column_thermodynamic_diagnostics!(thermodynamics)

    return nothing
end

