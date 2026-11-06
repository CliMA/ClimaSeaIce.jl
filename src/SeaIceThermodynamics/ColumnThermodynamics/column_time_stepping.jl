# The bottom flux is cached so that a stateful external flux is evaluated once per step.
@kernel function _cache_column_bottom_external_flux!(Qb, fields, grid, Qe, clock, model_fields)
    i, j = @index(Global, NTuple)
    @inbounds Qb[i, j, 1] = getflux(Qe, i, j, grid, fields.temperature[i, j, 1], clock, model_fields)
end

function cached_external_heat_fluxes(th::ColumnEnergyThermodynamics, Qe, clock, model_fields)
    grid = th.fields.internal_energy.grid
    cache = th.auxiliary.bottom_external_flux
    launch!(architecture(grid), grid, :xy, _cache_column_bottom_external_flux!, cache, th.fields, grid, Qe.bottom, clock, model_fields)
    return (top = Qe.top, bottom = cache)
end

function column_surface_energy_solve!(top_boundary, th::ColumnEnergyThermodynamics, Qe, clock, fields, hc, Δt)
    assemble_column_energy_system!(th, Qe, clock, fields, hc, Δt)
    solve_column_energy_system!(th)
    return nothing
end

# Every outer iteration restarts from the start-of-step enthalpy, as in Icepack `temperature_changes`.
function column_surface_energy_solve!(::MeltingConstrainedFluxBalance, thermodynamics::ColumnEnergyThermodynamics,
                                      external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt; max_iterations = 20)
    grid = thermodynamics.fields.internal_energy.grid
    tolerance = sqrt(eps(eltype(grid)))
    Tₛ = interior(thermodynamics.auxiliary.surface_temperature)
    E = interior(thermodynamics.fields.internal_energy)
    Nz = size(grid, 3)

    E⁰ = copy(E)
    interior(thermodynamics.auxiliary.surface_start_energy) .= view(E⁰, :, :, Nz:Nz)
    Tₛ⁻ = similar(Tₛ)

    for iteration in 1:max_iterations
        Tₛ⁻ .= Tₛ
        E .= E⁰
        compute_column_thermodynamic_diagnostics!(thermodynamics)
        assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
        solve_column_energy_system!(thermodynamics)
        compute_column_thermodynamic_diagnostics!(thermodynamics)
        launch!(architecture(grid), grid, :xy, _update_column_surface_temperature!, thermodynamics.auxiliary, thermodynamics.fields, grid,
                external_heat_fluxes.top, thermodynamics.relation, clock, model_fields)
        maximum(abs, Tₛ .- Tₛ⁻) <= tolerance && break
    end

    return nothing
end

function column_energy_time_step!(thermodynamics::ColumnEnergyThermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_transport_coefficients!(thermodynamics)
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
