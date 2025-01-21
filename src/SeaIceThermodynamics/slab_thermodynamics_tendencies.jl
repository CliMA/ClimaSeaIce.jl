using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: bottom_temperature, top_surface_temperature
import ClimaSeaIce.SeaIceThermodynamics: thickness_thermodynamic_tendency

@inline function thickness_thermodynamic_tendency(i, j, k, grid,
                                                  ice_thickness,
                                                  ice_concentration,
                                                  ice_consolidation_thickness,
                                                  thermodynamics::SlabSeaIceThermodynamics,
                                                  top_external_heat_flux,
                                                  bottom_external_heat_flux,
                                                  clock, model_fields)

    phase_transitions = thermodynamics.phase_transitions

    top_heat_bc = thermodynamics.heat_boundary_conditions.top
    bottom_heat_bc = thermodynamics.heat_boundary_conditions.bottom
    liquidus = phase_transitions.liquidus

    Qi = thermodynamics.internal_heat_flux
    Qu = top_external_heat_flux
    Qb = bottom_external_heat_flux
    Tu = thermodynamics.top_surface_temperature

    @inbounds begin
        hᶜ = ice_consolidation_thickness[i, j, k]
        hᵢ = ice_thickness[i, j, k]
    end

    # Consolidation criteria
    consolidated_ice = hᵢ >= hᶜ

    # Determine top surface temperature. 
    # Does this really fit here?
    # This is updating the temperature inside the thermodynamics module
    if !isa(top_heat_bc, PrescribedTemperature) # update surface temperature?
        if consolidated_ice # slab is consolidated and has an independent surface temperature
            Tu⁻ = @inbounds Tu[i, j, k]
            Tuⁿ = top_surface_temperature(i, j, grid, top_heat_bc, Tu⁻, Qi, Qu, clock, model_fields)
        else # slab is unconsolidated and does not have an independent surface temperature
            Tuⁿ = bottom_temperature(i, j, grid, bottom_heat_bc, liquidus)
        end

        @inbounds Tu[i, j, k] = Tuⁿ
    end

    @inbounds Tuᵢ = Tu[i, j, k]

    Tbᵢ = bottom_temperature(i, j, grid, bottom_heat_bc, liquidus)
    ℰb = latent_heat(phase_transitions, Tbᵢ)
    ℰu = latent_heat(phase_transitions, Tuᵢ)

    # Retrieve fluxes
    Quᵢ = getflux(Qu, i, j, grid, Tuᵢ, clock, model_fields)
    Qiᵢ = getflux(Qi, i, j, grid, Tuᵢ, clock, model_fields)
    Qbᵢ = getflux(Qb, i, j, grid, Tuᵢ, clock, model_fields)

    # If ice is consolidated, compute tendency for an ice slab; otherwise
    # just add ocean fluxes from frazil ice formation or melting
    slushy_Gh = - Qbᵢ / ℰb

    # Upper (top) and bottom interface velocities
    wu = (Quᵢ - Qiᵢ) / ℰu # < 0 => melting
    wb = (Qiᵢ - Qbᵢ) / ℰb # < 0 => freezing

    slabby_Gh = wu + wb

    return ifelse(consolidated_ice, slabby_Gh, slushy_Gh)
end

