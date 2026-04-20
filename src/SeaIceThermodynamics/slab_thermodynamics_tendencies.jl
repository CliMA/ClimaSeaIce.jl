using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: bottom_temperature, top_surface_temperature
using Oceananigans

#####
##### Ice interior conductive flux
#####

# Given the wrapped `FluxFunction` for the ice layer, evaluate the
# conductive flux at the ice-top temperature. Both the bare-ice path and
# the layered path build the wrapper inline before calling here.
@inline function ice_interior_heat_flux(Qi_function,
                                        i, j, k, grid, Tui,
                                        consolidated_ice, clock, model_fields)
    # If the ice is consolidated, we use the internal heat flux.
    # Slab is unconsolidated, there is no internal heat flux (Qi -> ∞)
    return ifelse(consolidated_ice,
                  getflux(Qi_function, i, j, grid, Tui, clock, model_fields),
                  zero(grid))
end

#####
##### Pure ice melt/freeze tendency
#####

# Given an externally-determined ice-top temperature `Tui` and top external
# flux `top_effective_heat_flux` (which may already include snow-surface
# absorption), compute the volume tendency of the ice slab. No surface solve
# is performed here.
@inline function ice_melt_freeze_tendency(i, j, k, grid,
                                          ice_thermodynamics::SlabThermodynamics,
                                          phase_transitions,
                                          sea_ice_density,
                                          Qi_function,
                                          Tui,
                                          ice_thickness,
                                          ice_consolidation_thickness,
                                          top_effective_heat_flux,
                                          bottom_external_heat_flux,
                                          clock, model_fields)

    bottom_heat_bc = ice_thermodynamics.heat_boundary_conditions.bottom
    liquidus = phase_transitions.liquidus

    @inbounds hi = ice_thickness[i, j, k]
    @inbounds hc = ice_consolidation_thickness[i, j, k]
    @inbounds ρi = sea_ice_density[i, j, 1]

    consolidated_ice = hi ≥ hc

    Tbi = bottom_temperature(i, j, grid, bottom_heat_bc, liquidus)

    # Energy per unit volume of sea-ice: bulk-density × per-mass latent heat
    ℰb = ρi * latent_heat(phase_transitions, Tbi)
    ℰu = ρi * latent_heat(phase_transitions, Tui)

    # Retrieve fluxes
    Qui = getflux(top_effective_heat_flux, i, j, grid, Tui, clock, model_fields)
    Qbi = getflux(bottom_external_heat_flux, i, j, grid, Tui, clock, model_fields)
    Qii = ice_interior_heat_flux(Qi_function, i, j, k, grid, Tui,
                                 consolidated_ice, clock, model_fields)

    # Upper (top) and bottom interface velocities
    # wu < 0 => top melting (volume loss from top)
    # wb > 0 => bottom freezing (volume gain at bottom)
    wu = (Qui - Qii) / ℰu
    wb = (Qii - Qbi) / ℰb

    return wu + wb
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

    @inbounds hi = ice_thickness[i, j, k]
    @inbounds hc = ice_consolidation_thickness[i, j, k]
    @inbounds Si = model_fields.S[i, j, k]

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

    return ice_melt_freeze_tendency(i, j, k, grid,
                                    ice_thermodynamics,
                                    phase_transitions,
                                    sea_ice_density,
                                    Qi_function,
                                    Tui,
                                    ice_thickness, ice_consolidation_thickness,
                                    Qu, bottom_external_heat_flux,
                                    clock, model_fields)
end
