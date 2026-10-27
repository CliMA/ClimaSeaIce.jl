# The bottom flux is cached so that a stateful external flux is evaluated once per step.
@kernel function _cache_column_bottom_external_flux!(Qb, fields, grid, Qe, clock, model_fields)
    i, j = @index(Global, NTuple)
    @inbounds Qb[i, j, 1] = getflux(Qe, i, j, grid, fields.temperature[i, j, 1], clock, model_fields)
end

function cached_external_heat_fluxes(th::ColumnEnergyThermodynamics, Qe, clock, model_fields)
    grid = th.fields.internal_energy.grid
    cache = th.auxiliary.bottom_external_flux
    launch!(architecture(grid), grid, :xy, _cache_column_bottom_external_flux!, cache, th.fields, grid, Qe.bottom, clock, model_fields)
    return (top = Qe.top, bottom = cached_column_bottom_flux(Qe.bottom, cache))
end

@inline cached_column_bottom_flux(bottom, cache) = cache
@inline cached_column_bottom_flux(bottom::ColumnBasalHeatFlux, cache) =
    ColumnBasalHeatFlux(cache, bottom.interface_cell, bottom.frazil)

@inline update_column_surface_temperature!(top_boundary, thermodynamics, external_flux, clock, model_fields) = nothing

function update_column_surface_temperature!(::MeltingConstrainedFluxBalance, thermodynamics, external_flux, clock, model_fields)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xy, _update_column_surface_temperature!, thermodynamics.auxiliary,
            thermodynamics.fields, grid, external_flux, thermodynamics.relation, clock, model_fields)
    return nothing
end

# A 90-minute solve can alternate between a cold surface and the melting cap:
# each branch is implicit, but selecting the branch from the previous candidate
# is an active-set iteration. Relax only columns whose branch changes.
@kernel function _relax_column_energy_at_surface_switch!(E, E_previous, Tₛ, Tₛ_previous, fields, relation, grid)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)
    @inbounds begin
        Tₘ = melting_temperature(relation.phase_transitions.liquidus, fields.bulk_salinity[i, j, Nz])
        was_melting = Tₛ_previous[i, j, 1] >= Tₘ
        is_melting = Tₛ[i, j, 1] >= Tₘ
        if was_melting != is_melting
            E[i, j, k] = (E_previous[i, j, k] + E[i, j, k]) / 2
        end
    end
end

# The enthalpy-to-temperature relation is nonlinear near the melting point.
# A Picard update can overshoot even when the conductive matrix itself is
# implicit. Dampen a layer only when its last two updates reverse direction;
# this leaves convergent and linear solves untouched.
@kernel function _relax_oscillating_column_energy!(E, E_previous, E_before_previous)
    i, j, k = @index(Global, NTuple)
    @inbounds begin
        previous_update = E_previous[i, j, k] - E_before_previous[i, j, k]
        current_update = E[i, j, k] - E_previous[i, j, k]
        if previous_update * current_update < 0
            E[i, j, k] = (E[i, j, k] + E_previous[i, j, k]) / 2
        end
    end
end

# Each solve uses the same start-of-step enthalpy on the right-hand side. Its
# matrix and affine fluxes are linearized around the previous candidate, so
# both interior conduction and the surface balance converge at the new time.
function column_surface_energy_solve!(top_boundary, thermodynamics::ColumnEnergyThermodynamics,
                                      external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt; max_iterations = 64)
    grid = thermodynamics.fields.internal_energy.grid
    # The skin root is solved to an absolute flux residual of 1e-6 W/m².
    # Its corresponding temperature can jitter by O(1e-7) K after the column
    # energy has converged. Do not mistake that root precision for divergence.
    surface_tolerance = max(sqrt(eps(eltype(grid))), 1e-6)
    Tₛ = interior(thermodynamics.auxiliary.surface_temperature)
    E = interior(thermodynamics.fields.internal_energy)
    Nz = size(grid, 3)

    E⁰ = copy(E)
    interior(thermodynamics.auxiliary.surface_start_energy) .= view(E⁰, :, :, Nz:Nz)
    E⁻ = similar(E)
    E⁻⁻ = similar(E)
    Tₛ⁻ = copy(Tₛ)
    energy_change = Inf
    surface_change = Inf

    for iteration in 1:max_iterations
        E⁻ .= E
        Tₛ⁻ .= Tₛ
        compute_column_thermodynamic_diagnostics!(thermodynamics)
        compute_column_transport_coefficients!(thermodynamics)
        assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt; initial_energy = E⁰)
        solve_column_energy_system!(thermodynamics)
        iteration > 1 && launch!(architecture(grid), grid, :xyz, _relax_oscillating_column_energy!, E, E⁻, E⁻⁻)
        compute_column_thermodynamic_diagnostics!(thermodynamics)
        update_column_surface_temperature!(top_boundary, thermodynamics, external_heat_fluxes.top, clock, model_fields)
        if top_boundary isa MeltingConstrainedFluxBalance
            launch!(architecture(grid), grid, :xyz, _relax_column_energy_at_surface_switch!,
                    E, E⁻, Tₛ, Tₛ⁻, thermodynamics.fields, thermodynamics.relation, grid)
            compute_column_thermodynamic_diagnostics!(thermodynamics)
            update_column_surface_temperature!(top_boundary, thermodynamics, external_heat_fluxes.top, clock, model_fields)
        end
        energy_tolerance = max(1e-8, 4eps(eltype(grid))) * max(one(eltype(grid)), maximum(abs, E))
        energy_change = maximum(abs, E .- E⁻)
        surface_change = maximum(abs, Tₛ .- Tₛ⁻)
        energy_converged = energy_change <= energy_tolerance
        surface_converged = surface_change <= surface_tolerance
        energy_converged && surface_converged && return nothing
        E⁻⁻ .= E⁻
    end

    ΔE = Array(E .- E⁻)
    eidx = argmax(abs.(ΔE))
    ΔTs = Array(Tₛ .- Tₛ⁻)
    sidx = argmax(abs.(ΔTs))
    T = Array(interior(thermodynamics.fields.temperature))
    TE = Array(interior(thermodynamics.fields.temperature_energy_derivative))
    S = Array(interior(thermodynamics.fields.bulk_salinity))
    Efinal = Array(E)
    error("Column energy/surface solve did not converge within $max_iterations iterations " *
          "(energy change $energy_change, tolerance $energy_tolerance at $eidx: E=$(Efinal[eidx]), T=$(T[eidx]), " *
          "dT/dE=$(TE[eidx]), S=$(S[eidx]); " *
          "surface change $surface_change, tolerance $surface_tolerance at $sidx: Ts=$(Array(Tₛ)[sidx]))")
end

function column_energy_time_step!(thermodynamics::ColumnEnergyThermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_shortwave_flux!(thermodynamics)
    column_surface_energy_solve!(thermodynamics.heat_boundary_conditions.top, thermodynamics, external_heat_fluxes, clock, model_fields,
                                 consolidation_thickness, Δt)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    column_salinity_time_step!(thermodynamics, Δt)
    return nothing
end

# The column carries its own enthalpy and interface heights, which RK substeps do not reset, so it steps once per full step.
@inline apply_column_thermodynamics_this_stage(timestepper, clock) = true
@inline apply_column_thermodynamics_this_stage(timestepper::SplitRungeKuttaTimeStepper, clock) = clock.stage == timestepper.Nstages

function thermodynamic_time_step!(model, th::ColumnEnergyThermodynamics, ::Nothing, Δt)
    apply_column_thermodynamics_this_stage(model.timestepper, model.clock) || return nothing
    Qe = cached_external_heat_fluxes(th, model.external_heat_fluxes, model.clock, fields(model))
    column_energy_time_step!(th, Qe, model.clock, fields(model), model.ice_consolidation_thickness, Δt)
    column_stefan_volume_update!(th, model, Qe, Δt)
    return nothing
end
