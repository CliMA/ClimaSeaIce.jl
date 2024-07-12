using Oceananigans.Advection
using ClimaSeaIce.SeaIceThermodynamics: ice_thickness_growth

@kernel function _compute_tracer_tendencies!(Gⁿ, ice_thickness,
                                             grid,
                                             clock,
                                             velocities,
                                             advection,
                                             concentration,
                                             thermodynamics,
                                             top_external_heat_flux,
                                             bottom_external_heat_flux,
                                             h_forcing,
                                             model_fields)

    i, j, k = @index(Global, NTuple)
    
    @inbounds Gⁿ.h[i, j, k] = ice_thickness_tendency(i, j, k, grid, clock,
                                                     velocities,
                                                     advection,
                                                     ice_thickness,
                                                     concentration,
                                                     thermodynamics,
                                                     top_external_heat_flux,
                                                     bottom_external_heat_flux,
                                                     h_forcing,
                                                     model_fields)
     
    @inbounds Gⁿ.ℵ[i, j, k] = - horizontal_div_Uc(i, j, k, grid, advection, velocities, concentration)
end

# Thickness change due to accretion and melting, restricted by minimum allowable value
function ice_thickness_tendency(i, j, k, grid, clock,
                                velocities,
                                advection,
                                ice_thickness,
                                concentration,
                                thermodynamics,
                                top_external_heat_flux,
                                bottom_external_heat_flux,
                                h_forcing,
                                model_fields)

    Gh_advection = - div_Uℵh(i, j, k, grid, advection, velocities, concentration, ice_thickness)

    Gh_growth = ice_thickness_growth(i, j, k, grid, 
                                     ice_thickness, 
                                     concentration,
                                     thermodynamics,
                                     top_external_heat_flux,
                                     bottom_external_heat_flux,
                                     clock, model_fields)

    
    # Compute forcing
    Fh = zero(grid) #h_forcing(i, j, grid, clock, model_fields)

    return Gh_advection + Gh_growth + Fh
end
