using Oceananigans.Advection
using ClimaSeaIce.SeaIceThermodynamics: thickness_thermodynamic_tendency

# Thickness change due to accretion and melting, restricted by minimum allowable value
function ice_thickness_tendency(i, j, k, grid, clock,
                                velocities,
                                advection,
                                ice_thickness,
                                ice_concentration,
                                thermodynamics,
                                top_external_heat_flux,
                                bottom_external_heat_flux,
                                h_forcing,
                                model_fields)

    Gh_advection = - div_Uc(i, j, k, grid, advection, velocities, ice_thickness)

    Gh_growth = thickness_thermodynamic_tendency(i, j, k, grid, 
                                                 ice_thickness, 
                                                 ice_concentration,
                                                 thermodynamics,
                                                 top_external_heat_flux,
                                                 bottom_external_heat_flux,
                                                 clock, model_fields)

    
    # Compute forcing
    Fh = zero(grid) #h_forcing(i, j, grid, clock, model_fields)

    return Gh_advection + Gh_growth + Fh
end
