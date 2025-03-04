using Oceananigans.Advection
using ClimaSeaIce.SeaIceThermodynamics: thickness_thermodynamic_tendency
using ClimaSeaIce.SeaIceMomentumEquations: compute_momentum_tendencies!

function compute_tendencies!(model::SIM, Δt)
    compute_tracer_tendencies!(model)
    compute_momentum_tendencies!(model, model.dynamics, Δt)
    return nothing
end

function compute_tracer_tendencies!(model::SIM)
    grid = model.grid
    arch = architecture(grid)
   
    launch!(arch, grid, :xy,
            _compute_tracer_tendencies!,
            model.timestepper.Gⁿ,
            model.ice_thickness,
            grid,
            model.clock,
            model.velocities,
            model.advection,
            model.ice_concentration,
            model.ice_consolidation_thickness,
            model.ice_thermodynamics,
            model.external_heat_fluxes.top,
            model.external_heat_fluxes.bottom,
            model.forcing.h,
            fields(model))

    return nothing
end

@kernel function _compute_tracer_tendencies!(Gⁿ, ice_thickness,
                                             grid,
                                             clock,
                                             velocities,
                                             advection,
                                             ice_concentration,
                                             ice_consolidation_thickness,
                                             thermodynamics,
                                             top_external_heat_flux,
                                             bottom_external_heat_flux,
                                             h_forcing,
                                             model_fields)

    i, j = @index(Global, NTuple)
    
    @inbounds Gⁿ.h[i, j, 1] = ice_thickness_tendency(i, j, 1, grid, clock,
                                                     velocities,
                                                     advection,
                                                     ice_thickness,
                                                     ice_concentration,
                                                     ice_consolidation_thickness,
                                                     thermodynamics,
                                                     top_external_heat_flux,
                                                     bottom_external_heat_flux,
                                                     h_forcing,
                                                     model_fields)
     
    @inbounds Gⁿ.ℵ[i, j, 1] = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, ice_concentration)
end

# Thickness change due to accretion and melting, restricted by minimum allowable value
function ice_thickness_tendency(i, j, k, grid, clock,
                                velocities,
                                advection,
                                ice_thickness,
                                ice_concentration,
                                ice_consolidation_thickness,
                                thermodynamics,
                                top_external_heat_flux,
                                bottom_external_heat_flux,
                                h_forcing,
                                model_fields)

    Gh_advection = - horizontal_div_Uc(i, j, k, grid, advection, velocities, ice_thickness) # div_Uℵh(i, j, k, grid, advection, velocities, ice_concentration, ice_thickness)

    Gh_thermodynamics = thickness_thermodynamic_tendency(i, j, k, grid, 
                                                         ice_thickness, 
                                                         ice_concentration,
                                                         ice_consolidation_thickness,
                                                         thermodynamics,
                                                         top_external_heat_flux,
                                                         bottom_external_heat_flux,
                                                         clock, model_fields)

    
    # Compute forcing
    Fh = zero(grid) #h_forcing(i, j, grid, clock, model_fields)

    return Gh_advection + Gh_thermodynamics + Fh 
end

