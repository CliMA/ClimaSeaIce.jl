using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: bottom_temperature, top_surface_temperature
using Oceananigans

@inline function thermodynamic_tendency(i, j, k, grid,
                                        ice_thermodynamics::SlabThermodynamics,
                                        ice_thickness,
                                        ice_concentration,
                                        ice_consolidation_thickness,
                                        top_external_heat_flux,
                                        bottom_external_heat_flux,
                                        clock, model_fields)

    phase_transitions = ice_thermodynamics.phase_transitions

    top_heat_bc = ice_thermodynamics.heat_boundary_conditions.top
    bottom_heat_bc = ice_thermodynamics.heat_boundary_conditions.bottom
    liquidus = phase_transitions.liquidus

    Qi = ice_thermodynamics.internal_heat_flux
    Qu = top_external_heat_flux
    Qb = bottom_external_heat_flux
    Tu = ice_thermodynamics.top_surface_temperature

    @inbounds begin
        hi = ice_thickness[i, j, k]
        hc = ice_consolidation_thickness[i, j, k]
        Si = model_fields.S[i, j, k]
    end

    @inbounds Tui = Tu[i, j, k]

    consolidated_ice = hi ≥ hc

    # Determine top surface temperature. 
    # Does this really fit here?
    # This is updating the temperature inside the ice_thermodynamics module
    if !isa(top_heat_bc, PrescribedTemperature) # update surface temperature?
        if consolidated_ice # slab is consolidated and has an independent surface temperature
            Tu⁻ = @inbounds Tu[i, j, k]
            Tuⁿ = top_surface_temperature(i, j, grid, top_heat_bc, Tu⁻, Qi, Qu, clock, model_fields)
            # We cap by melting temperature
            Tuₘ = melting_temperature(liquidus, Si)
            Tuⁿ = min(Tuⁿ, Tuₘ)
        else # slab is unconsolidated and does not have an independent surface temperature
            Tuⁿ = bottom_temperature(i, j, grid, bottom_heat_bc, liquidus)
        end

        @inbounds Tu[i, j, k] = Tuⁿ
    end

    @inbounds Tui = Tu[i, j, k]

    Tbi = bottom_temperature(i, j, grid, bottom_heat_bc, liquidus)
    ℰb = latent_heat(phase_transitions, Tbi)
    ℰu = latent_heat(phase_transitions, Tui)

    # Retrieve fluxes
    Qui = getflux(Qu, i, j, grid, Tui, clock, model_fields)
    Qbi = getflux(Qb, i, j, grid, Tui, clock, model_fields)

    if consolidated_ice # If the ice is consolidated, we use the internal heat flux
        Qii = getflux(Qi, i, j, grid, Tui, clock, model_fields)
    else # Slab is unconsolidated, there is no internal heat flux (Qi -> ∞)
        Qii = zero(grid)
    end
    
    # Upper (top) and bottom interface velocities
    # wu < 0 => top melting (volume loss from top)
    # wb > 0 => bottom freezing (volume gain at bottom)
    wu = (Qui - Qii) / ℰu
    wb = (Qii - Qbi) / ℰb

    return wu + wb
end
