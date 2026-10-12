using .HeatBoundaryConditions: bottom_temperature, top_surface_temperature

#####
##### Pure ice melt/freeze tendency
#####

# Given an externally-determined ice-top temperature `Tui` and top external
# flux `top_effective_heat_flux` (which may already include snow-surface
# absorption), compute the volume tendency of the ice slab. No surface solve
# is performed here.
@inline function ice_melt_freeze_tendency(i, j, k, grid,
                                          phase_transitions,
                                          sea_ice_density,
                                          Tui,
                                          top_effective_heat_flux,
                                          bottom_external_heat_flux,
                                          clock, model_fields)

    @inbounds ρi = sea_ice_density[i, j, 1]

    # The slab stores ρi ℒ₀ per unit volume and the internal flux only moves energy between its two interfaces
    ℰ = ρi * phase_transitions.reference_latent_heat

    Qui = getflux(top_effective_heat_flux, i, j, grid, Tui, clock, model_fields)
    Qbi = getflux(bottom_external_heat_flux, i, j, grid, Tui, clock, model_fields)

    return (Qui - Qbi) / ℰ
end

#####
##### Top-level tendency with surface solve (bare-ice entry point)
#####

@inline function thermodynamic_tendency(i, j, k, grid,
                                        ice_thermodynamics::SlabThermodynamics,
                                        phase_transitions,
                                        sea_ice_density,
                                        ice_thickness,
                                        ice_concentration,
                                        ice_consolidation_thickness,
                                        top_external_heat_flux,
                                        bottom_external_heat_flux,
                                        clock, model_fields)

    top_heat_bc = ice_thermodynamics.heat_boundary_conditions.top
    bottom_heat_bc = ice_thermodynamics.heat_boundary_conditions.bottom
    liquidus = phase_transitions.liquidus

    # Build the internal-flux wrapper inline using the model's shared liquidus.
    Qi_function = internal_flux_function(ice_thermodynamics.internal_heat_flux,
                                         liquidus, bottom_heat_bc)
    Qu = top_external_heat_flux
    Tu = ice_thermodynamics.top_surface_temperature

    @inbounds begin
        hi = ice_thickness[i, j, k]
        hc = ice_consolidation_thickness[i, j, k]
        Si = model_fields.S[i, j, k]
    end

    consolidated_ice = hi ≥ hc

    # Determine top surface temperature.
    # Does this really fit here?
    # This is updating the temperature inside the ice_thermodynamics module
    if !isa(top_heat_bc, PrescribedTemperature) # update surface temperature?
        if consolidated_ice # slab is consolidated and has an independent surface temperature
            @inbounds Tu⁻ = Tu[i, j, k]
            Tuⁿ = top_surface_temperature(i, j, grid, top_heat_bc, Tu⁻, Qi_function, Qu, clock, model_fields)
            # We cap by melting temperature
            Tuₘ = melting_temperature(liquidus, Si)
            Tuⁿ = min(Tuⁿ, Tuₘ)
        else # slab is unconsolidated and does not have an independent surface temperature
            Tuⁿ = bottom_temperature(i, j, grid, bottom_heat_bc, liquidus)
        end
        @inbounds Tu[i, j, k] = Tuⁿ
    end

    @inbounds Tui = Tu[i, j, k]

    # Evaluate the external flux closures exactly once at the converged Tui and pass the
    # scalar values down to `ice_melt_freeze_tendency` (which accepts Number via `getflux`).
    Qui = getflux(Qu, i, j, grid, Tui, clock, model_fields)
    Qbi = getflux(bottom_external_heat_flux, i, j, grid, Tui, clock, model_fields)

    return ice_melt_freeze_tendency(i, j, k, grid,
                                    phase_transitions,
                                    sea_ice_density,
                                    Tui,
                                    Qui, Qbi,
                                    clock, model_fields)
end
