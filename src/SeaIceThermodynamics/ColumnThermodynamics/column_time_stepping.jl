# Evaluate the bottom external flux exactly once per step into a field, so a stateful coupler (e.g. an evolving
# lake/ocean) fires its side effect only once even though the energy solve and the volume update both read it.
# `getflux` of the resulting 2D field is a plain array read. The top flux stays live because a melting surface
# solve must evaluate it at trial temperatures, exactly as SlabThermodynamics's root-finder does.
@kernel function _cache_column_bottom_external_flux!(Qb, fields, grid, Qe, clock, model_fields)
    i, j = @index(Global, NTuple)
    T = @inbounds fields.temperature[i, j, 1]
    @inbounds Qb[i, j, 1] = getflux(Qe, i, j, grid, T, clock, model_fields)
end

function cached_external_heat_fluxes(th::ColumnEnergyThermodynamics, Qe, clock, model_fields)
    grid = th.fields.internal_energy.grid
    cache = th.auxiliary.bottom_external_flux
    launch!(architecture(grid), grid, :xy, _cache_column_bottom_external_flux!, cache, th.fields, grid, Qe.bottom, clock, model_fields)
    return (top = Qe.top, bottom = cache)
end

# The surface energy solve is a single tridiagonal pass for direct-flux and Dirichlet tops; a
# `MeltingConstrainedFluxBalance` top adds the outer surface-temperature iteration around it.
function column_surface_energy_solve!(th::ColumnEnergyThermodynamics, Qe, clock, fields, hc, Δt)
    column_surface_energy_solve!(th.heat_boundary_conditions.top, th, Qe, clock, fields, hc, Δt)
    return nothing
end

function column_surface_energy_solve!(top_boundary, th::ColumnEnergyThermodynamics, Qe, clock, fields, hc, Δt)
    assemble_column_energy_system!(th, Qe, clock, fields, hc, Δt)
    solve_column_energy_system!(th)
    return nothing
end

function column_surface_energy_solve!(top_boundary::MeltingConstrainedFluxBalance,
                                      thermodynamics::ColumnEnergyThermodynamics,
                                      external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt;
                                      max_iterations = 20)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)
    tolerance = sqrt(eps(eltype(grid)))
    surface_temperature = thermodynamics.auxiliary.surface_temperature
    internal_energy = thermodynamics.fields.internal_energy

    # Each outer pass re-solves from the same start-of-step enthalpy with the latest surface coupling, exactly as
    # Icepack `temperature_changes` rebuilds the system from the initial temperature every iteration. Restoring it
    # (rather than accumulating onto the previous solve) is what keeps the lagged top flux applied once per step.
    start_of_step_energy = copy(interior(internal_energy))

    # Stash the top-cell start-of-step enthalpy so the basal/surface Stefan residual in the volume update is computed
    # against the same enthalpy the melt cap used, not the post-solve enthalpy (which the cap already drove to melt).
    Nz = size(grid, 3)
    interior(thermodynamics.auxiliary.surface_start_energy)[:, :, 1] .= @view start_of_step_energy[:, :, Nz]

    for iteration in 1:max_iterations
        previous = Array(interior(surface_temperature))
        interior(internal_energy) .= start_of_step_energy
        compute_column_thermodynamic_diagnostics!(thermodynamics)
        assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
        solve_column_energy_system!(thermodynamics)
        compute_column_thermodynamic_diagnostics!(thermodynamics)
        launch!(arch, grid, :xy, _update_column_surface_temperature!,
                thermodynamics.auxiliary, thermodynamics.fields, grid,
                external_heat_fluxes.top, thermodynamics.relation, clock, model_fields)
        maximum(abs.(Array(interior(surface_temperature)) .- previous)) <= tolerance && break
    end

    return nothing
end

# Standalone solves treat every column as consolidated (`nothing`); the coupled step passes the model's
# consolidation thickness so thin columns take the slab-balance regime.
column_energy_time_step!(thermodynamics::ColumnEnergyThermodynamics, external_heat_fluxes, clock, model_fields, Δt) =
    column_energy_time_step!(thermodynamics, external_heat_fluxes, clock, model_fields, nothing, Δt)

function column_energy_time_step!(thermodynamics::ColumnEnergyThermodynamics,
                                  external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_transport_coefficients!(thermodynamics)
    # The shortwave deposition depends only on grid geometry and the absorption closure, so it is computed once per
    # step here rather than inside every (possibly iterated) energy-system assembly.
    compute_column_shortwave_flux!(thermodynamics)
    column_surface_energy_solve!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    column_salinity_time_step!(thermodynamics, Δt)

    return nothing
end

# The column carries its own prognostic enthalpy and moving vertical metric (unlike the slab, which derives
# everything from the thickness). The RK substeps reset the thickness from the start-of-step cache but not the
# column's enthalpy/grid, so applying the column physics on every substep would re-advance that carried state from
# each substep's partially-stepped values (e.g. re-melting a surface cell already driven to complete melt). Apply it
# once per full step, on the final RK stage where `Δτ = Δt` and the earlier stages leave the enthalpy/grid intact.
@inline apply_column_thermodynamics_this_stage(timestepper, clock) = true
@inline apply_column_thermodynamics_this_stage(timestepper::SplitRungeKuttaTimeStepper, clock) =
    clock.stage == timestepper.Nstages

function thermodynamic_time_step!(model, th::ColumnEnergyThermodynamics, ::Nothing, Δt)

    apply_column_thermodynamics_this_stage(model.timestepper, model.clock) || return nothing

    # Evaluate stateful external fluxes once per step, then drive both the energy solve and the volume update from
    # the cached values so their side effects (e.g. an evolving lake) fire exactly once — as in SlabThermodynamics.
    Qe = cached_external_heat_fluxes(th, model.external_heat_fluxes, model.clock, fields(model))
    column_energy_time_step!(th, Qe, model.clock, fields(model), model.ice_consolidation_thickness, Δt)
    column_stefan_volume_update!(th, model, Qe, Δt)
    return nothing
end

