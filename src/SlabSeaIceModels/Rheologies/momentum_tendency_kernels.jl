""" tendency computation for the ice velocities. Used only in case of
an explicit (or semi-implicit) rheology """
@kernel function _compute_momentum_tendencies!(tendencies,
                                               grid,
                                               clock,
                                               velocities,
                                               ocean_velocities,
                                               coriolis,
                                               thickness,
                                               rheology,
                                               top_u_stress,
                                               top_v_stress,
                                               forcing,
                                               model_fields)

    i, j = @index(Global, NTuple)

    Gu = tendencies.u
    Gv = tendencies.v

    τua = top_u_stress
    τva = top_v_stress
    
    momentum_args = (i, j, grid, clock, velocities, ocean_velocities, coriolis, rheology)

    τuo, τvo = ice_ocean_stress(i, j, ocean_velocities, velocities, thickness)

    @inbounds Gu[i, j, 1] = u_velocity_tendency(momentum_args..., τua, nothing, model_fields) + τuo
    @inbounds Gv[i, j, 1] = v_velocity_tendency(momentum_args..., τva, nothing, model_fields) + τvo
end

""" tendency computation for the ice u-velocity """
function u_velocity_tendency(i, j, grid, clock,
                             velocities, 
                             ocean_velocities,
                             coriolis,
                             rheology,
                             u_top_stress,
                             u_forcing,
                             model_fields)

    u,  v  = velocities
    uₒ, vₒ = ocean_velocities

    relative_u = DifferenceOfArrays(u, uₒ)
    relative_v = DifferenceOfArrays(v, vₒ)

    @inbounds top_boundary = u_top_stress[i, j, 1]

    relative_velocities = (; u = relative_u, v = relative_v)
    
    return ( - x_f_cross_U(i, j, 1, grid, coriolis, relative_velocities)
             + boundary
             + x_stress_divergence(i, j, 1, rheology) )
end

""" tendency computation for the ice v-velocity """
function v_velocity_tendency(i, j, grid, clock,
                             velocities,
                             ocean_velocities, 
                             coriolis,
                             rheology,
                             v_top_stress,
                             v_forcing,
                             model_fields)

    u,  v  = velocities
    uₒ, vₒ = ocean_velocities

    relative_u = DifferenceOfArrays(u, uₒ)
    relative_v = DifferenceOfArrays(v, vₒ)

    relative_velocities = (; u = relative_u, v = relative_v)
    
    @inbounds top_boundary = v_top_stress[i, j, 1]

    return  ( - y_f_cross_U(i, j, 1, grid, coriolis, relative_velocities)
              + boundary
              + y_stress_divergence(i, j, 1, rheology) )
end
