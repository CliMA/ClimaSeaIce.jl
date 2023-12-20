
""" stepping the ice u-velocity using a forward leap-frog scheme """
@kernel function _u_velocity_step!(velocities, grid, Δt, 
                                   clock,
                                   ocean_velocities,
                                   coriolis,
                                   rheology,
                                   thickness,
                                   u_top_stress,
                                   u_forcing,
                                   model_fields)

    i, j = @index(Global, NTuple)

    uₒ, vₒ = ocean_velocities

    @inbounds τuₐ = u_top_stress[i, j, 1]
    τiₒ = ice_ocean_implicit_stress(i, j, ocean_velocities, velocities, thickness)

    @inbounds Gᵘ = ( - x_f_cross_U(i, j, 1, grid, coriolis, velocities) 
                     + τuₐ
                     - τiₒ * uₒ[i, j, 1]
                     + x_stress_divergence(i, j, 1, rheology) )

    @inbounds u[i, j, 1] += Δt * Gᵘ
    @inbounds u[i, j, 1] /= (1 + Δt * τiₒ)
end

""" tendency computation for the ice v-velocity """
@kernel function _v_velocity_step!(velocities, grid, Δt, 
                                   clock,
                                   ocean_velocities, 
                                   coriolis,
                                   rheology,
                                   v_top_stress,
                                   v_forcing,
                                   model_fields)

    i, j = @index(Global, NTuple)

    @inbounds τva = v_top_stress[i, j, 1]
    τiₒ = ice_ocean_implicit_stress(i, j, ocean_velocities, velocities, thickness)

    @inbounds Gᵛ = ( - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
                     + τva
                     - τiₒ * uₒ[i, j, 1]
                     + y_stress_divergence(i, j, 1, rheology) )

    @inbounds v[i, j, 1] += Δt * Gᵘ
    @inbounds v[i, j, 1] /= (1 + Δt * τiₒ)
end


# The stress is treated implicitly 
# i.e:
#
#      Cᴰρₒ
# τₒ = ---- || u - uₒ ||ⁿ * (uⁿ⁺¹ - uₒⁿ)
#      ρᵢ h
#
@inline function ice_ocean_implicit_stress(i, j, vel_oce, vel_ice, h)
    FT = eltype(h)

    # To put somewhere else??
    ρₒ = 1024
    ρᵢ = 917
    Cₒ = convert(FT, 1e-3)
    uᵢ, vᵢ = vel_ice
    uₒ, vₒ = vel_oce

    @inbounds begin
        δu    = uₒ[i, j, 1] - uᵢ[i, j, 1]
        δv    = vₒ[i, j, 1] - vᵢ[i, j, 1]
        δ     = sqrt(δu^2 + δv^2)
        τ_imp = ifelse(h[i, j, 1] > 0, Cₒ * ρₒ * δ / (h[i, j, 1] * ρᵢ), 0)
    end

    return τ_imp
end