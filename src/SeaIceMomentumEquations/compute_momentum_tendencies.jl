using Oceananigans.Utils

function compute_momentum_tendencies!(model, ::SeaIceMomentumEquations)
    
    ice_dynamics = model.ice_dynamics
    grid = model.grid

    args = (model.clock, 
            model.velocities, 
            model.ocean_free_surface, 
            ice_dynamics.coriolis, 
            ice_dynamics.rheology, 
            ice_dynamics.auxiliary_fields, 
            model.ice_thickness, 
            model.ice_concentration, 
            model.ice_density, 
            ice_dynamics.gravitational_acceleration)

    u_top_stress = model.external_momentum_stresses.top
    v_top_stress = model.external_momentum_stresses.top
    u_bottom_stress = model.external_momentum_stresses.bottom
    v_bottom_stress = model.external_momentum_stresses.bottom

    Gu = model.timestepper.Gⁿ.u
    Gv = model.timestepper.Gⁿ.v

    launch!(architecture(grid), grid, :xy, _compute_velocity_tendencies!, Gu, Gv, args,
            u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress)

    return nothing
end

@kernel function _compute_velocity_tendencies!(Gu, Gv, grid, args, u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress)
    i, j = @index(Global, NTuple)
    @inbounds Gu[i, j, k] = u_velocity_tendency(i, j, grid, args..., u_top_stress, u_bottom_stress)
    @inbounds Gv[i, j, k] = u_velocity_tendency(i, j, grid, args..., v_top_stress, v_bottom_stress)
end