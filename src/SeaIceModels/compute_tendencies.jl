using ClimaSeaIce.Advection: div_Uℵh, horizontal_div_Uc
using ClimaSeaIce.SeaIceThermodynamics: thickness_thermodynamic_tendency

function compute_tendencies!(model::SIM)
    compute_tracer_tendencies!(model)
    # compute_momentum_tendencies!(model)
    return nothing
end

function compute_tracer_tendencies!(model::SIM)
    grid = model.grid
    arch = architecture(grid)
   
    launch!(arch, grid, :xyz,
            _compute_tracer_tendencies!,
            model.timestepper.Gⁿ,
            model.ice_thickness,
            grid,
            model.clock,
            model.velocities,
            model.advection,
            model.ice_concentration,
            model.ice_thermodynamics,
            model.external_heat_fluxes.top,
            model.external_heat_fluxes.bottom,
            fields(model))

    return nothing
end

@kernel function _compute_tracer_tendencies!(Gⁿ, ice_thickness,
                                             grid,
                                             clock,
                                             velocities,
                                             advection,
                                             concentration,
                                             thermodynamics,
                                             top_external_heat_flux,
                                             bottom_external_heat_flux,
                                             model_fields)

    i, j, k = @index(Global, NTuple)

    Gh_advection = - div_Uℵh(i, j, k, grid, advection, velocities, ice_concentration, ice_thickness)

    Gh_thermodynamics = thickness_thermodynamic_tendency(i, j, k, grid, 
                                                         ice_thickness, 
                                                         ice_concentration,
                                                         thermodynamics,
                                                         top_external_heat_flux,
                                                         bottom_external_heat_flux,
                                                         clock, model_fields)
    
    @inbounds Gⁿ.h[i, j, k] = Gh_advection + Gh_thermodynamics
    @inbounds Gⁿ.ℵ[i, j, k] = - horizontal_div_Uc(i, j, k, grid, advection, velocities, concentration)
end

