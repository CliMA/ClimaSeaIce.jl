using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U
using ClimaSeaIce.SlabSeaIceModels.SlabSeaIceDynamics: Vᵢ

""" The beta coefficient for the leap-frog scheme """
@inline beta_coefficient(rheology, Δt) = 300

""" stepping the ice u-velocity using a forward leap-frog scheme """
@kernel function _u_velocity_step!(velocities, grid, Δt, 
                                   clock,
                                   ocean_velocities,
                                   coriolis,
                                   rheology,
                                   thickness,
                                   concentration,
                                   ice_density,
                                   u_top_stress,
                                   u_forcing,
                                   model_fields)

    i, j = @index(Global, NTuple)

    uᵢ, vᵢ = velocities
    uₒ, vₒ = ocean_velocities
    h  = thickness
    ℵ  = concentration
    ρᵢ = ice_density
    uⁿ = rheology.uⁿ

    # Ice mass interpolated on u points
    mᵢ = ℑxᶠᶜᶜ(i, j, 1, grid, Vᵢ, h, ℵ) * ρᵢ

    # relative ice velocities
    Δu = @inbounds uₒ[i, j, 1] - uᵢ[i, j, 1]
    Δv = ℑxyᶠᶜᶜ(i, j, 1, grid, vₒ) - ℑxyᶠᶜᶜ(i, j, 1, grid, vᵢ)

    # relative ice speed
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    β = beta_coefficient(rheology, Δt)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τuₐ = u_top_stress[i, j, 1]

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * uₒ - τₑₒ * uᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) 
    τₑₒ = 1e-3 * 1020 * Δ𝒰 / mᵢ

    @inbounds Gᵁ = ( - x_f_cross_U(i, j, 1, grid, coriolis, velocities) 
                     + τuₐ
                     + τₑₒ * uₒ[i, j, 1] # Explicit component of the ice-ocean stress
                     + x_internal_stress_divergence(i, j, grid, rheology) / mᵢ)

    # make sure we do not have NaNs!                 
    Gᵁ = ifelse(mᵢ > 0, Gᵁ, 0) 
    
    # Explicit step
    @inbounds uᵢ[i, j, 1] = (uᵢ[i, j, 1] * (β - 1) + Δt * Gᵁ + uⁿ[i, j, 1]) / β
    
    # Implicit component of the ice-ocean stress
    τᵢ = ifelse(mᵢ > 0, Δt * τₑₒ / β, 0)

    # Implicit step
    @inbounds uᵢ[i, j, 1] /= (1 + τᵢ) 
end

""" stepping the ice v-velocity using a forward leap-frog scheme """
@kernel function _v_velocity_step!(velocities, grid, Δt, 
                                   clock,
                                   ocean_velocities, 
                                   coriolis,
                                   rheology,
                                   thickness,
                                   concentration,
                                   ice_density,
                                   v_top_stress,
                                   v_forcing,
                                   model_fields)

    i, j = @index(Global, NTuple)

    uᵢ, vᵢ = velocities
    uₒ, vₒ = ocean_velocities
    h  = thickness
    ℵ  = concentration
    ρᵢ = ice_density
    vⁿ = rheology.vⁿ

    # Ice mass interpolated on v points
    mᵢ = ℑyᶜᶠᶜ(i, j, 1, grid, Vᵢ, h, ℵ) * ρᵢ
    
    # relative ice velocities
    Δu = ℑxyᶜᶠᶜ(i, j, 1, grid, uₒ) - ℑxyᶜᶠᶜ(i, j, 1, grid, uᵢ)
    Δv = @inbounds vₒ[i, j, 1] - vᵢ[i, j, 1]

    # relative ice speed
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    β = beta_coefficient(rheology, Δt)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τva = v_top_stress[i, j, 1] / mᵢ

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * vₒ - τₑₒ * vᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) / mᵢ
    τₑₒ = 1e-3 * 1020 * Δ𝒰 / mᵢ

    @inbounds Gⱽ = ( - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
                     + τva
                     + τₑₒ * vₒ[i, j, 1] # Explicit component of the ice-ocean stress
                     + y_internal_stress_divergence(i, j, grid, rheology) / mᵢ) 

    # make sure we do not have NaNs!
    Gⱽ = ifelse(mᵢ > 0, Gⱽ, 0) 

    # Explicit step
    @inbounds vᵢ[i, j, 1] = (vᵢ[i, j, 1] * (β - 1) + Δt * Gⱽ + vⁿ[i, j, 1]) / β

    # Implicit component of the ice-ocean stress
    τᵢ = ifelse(mᵢ > 0, Δt * τₑₒ / β, 0)

    # Implicit step
    @inbounds vᵢ[i, j, 1] /= (1 + τᵢ) 
end

# The ice-ocean stress is treated semi-implicitly 
# i.e:
#
#           Cᴰρₒ
# τₒ =    ------- || uₒ - uⁿ ||   * (uₒ - uⁿ⁺¹)
#         ρᵢ h ℵ
#      |-----------------------|
# τₒ =  τₑ₀ (explicit component)  * Δu   
#
#
