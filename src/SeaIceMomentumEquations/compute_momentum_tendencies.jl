using Oceananigans.Utils

# Compute the tendencies for the explicit momentum equations
function compute_momentum_tendencies!(model, ::ExplicitMomentumEquation)
    
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

    u_top_stress = model.external_momentum_stresses.top.u
    v_top_stress = model.external_momentum_stresses.top.v
    u_bottom_stress = model.external_momentum_stresses.bottom.u
    v_bottom_stress = model.external_momentum_stresses.bottom.v

    Gu = model.timestepper.Gⁿ.u
    Gv = model.timestepper.Gⁿ.v

    launch!(architecture(grid), grid, :xy, _compute_velocity_tendencies!, Gu, Gv, args,
            u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress)

    return nothing
end

@kernel function _compute_velocity_tendencies!(Gu, Gv, grid, args, u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress)
    i, j = @index(Global, NTuple)
    @inbounds Gu[i, j, k] = u_velocity_tendency(i, j, grid, args..., u_top_stress, u_bottom_stress)
    @inbounds Gv[i, j, k] = v_velocity_tendency(i, j, grid, args..., v_top_stress, v_bottom_stress)
end