using Oceananigans.Advection
using ClimaSeaIce.SeaIceThermodynamics: thickness_thermodynamic_tendency
using ClimaSeaIce.SeaIceMomentumEquations: compute_momentum_tendencies!

function compute_tendencies!(model::SIM)
    compute_tracer_tendencies!(model)
    compute_momentum_tendencies!(model, model.dynamics)

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
            model.ice_consolidation_thickness,
            model.ice_thermodynamics,
            model.tracers,
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
                                             tracers,
                                             top_external_heat_flux,
                                             bottom_external_heat_flux,
                                             h_forcing,
                                             model_fields)

    i, j, k = @index(Global, NTuple)
    
    @inbounds Gⁿ.h[i, j, k] = ice_thickness_tendency(i, j, k, grid, clock,
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
     
    @inbounds Gⁿ.ℵ[i, j, k] = - horizontal_div_Uc(i, j, k, grid, advection, velocities, ice_concentration)

    for n in eachindex(tracers) # all the other tracers just advect
        ϕ = @inbounds tracers[n]
        if !(ϕ isa ConstantField)
            @inbounds Gⁿ[n][i, j, k] = - horizontal_div_Uc(i, j, k, grid, advection, velocities, ϕ)
        end
    end
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
    Fh = h_forcing(i, j, grid, clock, model_fields)

    return Gh_advection + Gh_thermodynamics + Fh 
end

