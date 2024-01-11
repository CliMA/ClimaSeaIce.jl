using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U

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

    u, v   = velocities
    uₒ, vₒ = ocean_velocities

    @inbounds τuₐ = u_top_stress[i, j, 1]
    τiₒ = x_ice_ocean_implicit_stress(i, j, grid, ocean_velocities, velocities, thickness)

    @inbounds Gᵁ = ( - x_f_cross_U(i, j, 1, grid, coriolis, velocities) 
                     + τuₐ
                     + τiₒ * uₒ[i, j, 1]
                     + x_internal_stress_divergence(i, j, grid, rheology) )

    @inbounds u[i, j, 1] += Δt * Gᵁ
    @inbounds u[i, j, 1] /= (1 + Δt * τiₒ)
end

""" tendency computation for the ice v-velocity """
@kernel function _v_velocity_step!(velocities, grid, Δt, 
                                   clock,
                                   ocean_velocities, 
                                   coriolis,
                                   rheology,
                                   thickness,
                                   v_top_stress,
                                   v_forcing,
                                   model_fields)

    i, j = @index(Global, NTuple)

    u, v   = velocities
    uₒ, vₒ = ocean_velocities

    @inbounds τva = v_top_stress[i, j, 1]
    τiₒ = y_ice_ocean_implicit_stress(i, j, grid, ocean_velocities, velocities, thickness)

    @inbounds Gⱽ = ( - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
                     + τva
                     + τiₒ * vₒ[i, j, 1]
                     + y_internal_stress_divergence(i, j, grid, rheology) )

    @inbounds v[i, j, 1] += Δt * Gⱽ
    @inbounds v[i, j, 1] /= (1 + Δt * τiₒ)
end

# The ice-ocean stress is treated semi-implicitly 
# i.e:
#
#      Cᴰρₒ
# τₒ = ---- || uₒ - u ||ⁿ * (uₒⁿ - uⁿ⁺¹)
#      ρᵢ h
#     |------------------|
#   ice_ocean_implicit_stress
#
#
@inline function x_ice_ocean_implicit_stress(i, j, grid, vel_oce, vel_ice, h)
    FT = eltype(h)

    # To put somewhere else??
    ρₒ = 1024
    ρᵢ = 917
    Cₒ = convert(FT, 1e-3)

    uᵢ, vᵢ = vel_ice
    uₒ, vₒ = vel_oce

    @inbounds begin
        δu    = uₒ[i, j, 1] - uᵢ[i, j, 1]
        δv    = ℑxyᶠᶜᵃ(i, j, 1, grid, vₒ) - ℑxyᶠᶜᵃ(i, j, 1, grid, vᵢ)
        δ     = sqrt(δu^2 + δv^2)
        τ_imp = ifelse(h[i, j, 1] > 0, Cₒ * ρₒ * δ / (h[i, j, 1] * ρᵢ), 0)
    end

    return τ_imp
end

@inline function y_ice_ocean_implicit_stress(i, j, grid, vel_oce, vel_ice, h)
    FT = eltype(h)

    # To put somewhere else??
    ρₒ = 1024
    ρᵢ = 917
    Cₒ = convert(FT, 1e-3)

    uᵢ, vᵢ = vel_ice
    uₒ, vₒ = vel_oce

    @inbounds begin
        δu    = ℑxyᶜᶠᵃ(i, j, 1, grid, uₒ) - ℑxyᶜᶠᵃ(i, j, 1, grid, uᵢ)
        δv    = vₒ[i, j, 1] - vᵢ[i, j, 1]
        δ     = sqrt(δu^2 + δv^2)
        τ_imp = ifelse(h[i, j, 1] > 0, Cₒ * ρₒ * δ / (h[i, j, 1] * ρᵢ), 0)
    end

    return τ_imp
end