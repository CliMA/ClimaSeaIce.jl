using Oceananigans.Advection
using ClimaSeaIce.SeaIceThermodynamics: thickness_growth, bottom_ice_formation
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
    
    @inbounds hᵢ = ice_thickness[i, j, 1]
    @inbounds ℵᵢ = ice_concentration[i, j, 1]

    Gh⁺ = thickness_growth(i, j, 1, grid,
                           thermodynamics,
                           ice_thickness,
                           ice_concentration,
                           ice_consolidation_thickness,
                           top_external_heat_flux,
                           bottom_external_heat_flux,
                           clock, model_fields)

    GV⁺ = bottom_ice_formation(i, j, 1, grid, thermodynamics, bottom_external_heat_flux, clock, model_fields)

    Gℵ = GV⁺ / hᵢ
    Gh = Gh⁺ * ℵᵢ

    Guh = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, ice_thickness)
    Guℵ = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, ice_concentration)

    @inbounds Gⁿ.h[i, j, 1] = Gh + Guh
    @inbounds Gⁿ.ℵ[i, j, 1] = Gℵ + Guℵ
end

