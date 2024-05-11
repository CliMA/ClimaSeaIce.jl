using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U

@inline beta_coefficient(rheology, Δt) = rheology.substeps

""" stepping the ice u-velocity using a forward leap-frog scheme """
@kernel function _u_velocity_step!(velocities, grid, Δt, 
                                   clock,
                                   ocean_velocities,
                                   coriolis,
                                   rheology,
                                   thickness,
                                   concentration,
                                   u_top_stress,
                                   u_forcing,
                                   model_fields)

    i, j = @index(Global, NTuple)

    uᵢ, vᵢ = velocities
    uₒ, vₒ = ocean_velocities
    h  = thickness
    ℵ  = concentration
    ρᵢ = 917

    uⁿ = rheology.uⁿ

    mᵢ = @inbounds h[i, j, 1] * ℵ[i, j, 1] * ρᵢ
    Δu = @inbounds uₒ[i, j, 1] - uᵢ[i, j, 1]
    Δv = ℑxyᶠᶜᵃ(i, j, 1, grid, vₒ) - ℑxyᶠᶜᵃ(i, j, 1, grid, vᵢ)
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    β = beta_coefficient(rheology, Δt)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τuₐ = u_top_stress[i, j, 1]
    τₑₒ = 1e-3 * 1020 * Δ𝒰

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * uₒ - τₑₒ * uᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) 
    τₑₒ = 1e-3 * 1020 * Δ𝒰

    @inbounds Gᵁ = ( - x_f_cross_U(i, j, 1, grid, coriolis, velocities) 
                     + τuₐ
                     + τₑₒ * uₒ[i, j, 1] # Explicit component of the ice-ocean stress
                     + x_internal_stress_divergence(i, j, grid, rheology) ) / mᵢ

    # make sure we do not have NaNs!                 
    Gᵁ = ifelse(mᵢ > 0, Gᵁ, 0) 
    
    # Explicit step
    @inbounds uᵢ[i, j, 1] = (uᵢ[i, j, 1] * (β - 1) + Δt * Gᵁ + uⁿ[i, j, 1]) / β
    
    # Implicit component of the ice-ocean stress
    τᵢ = ifelse(mᵢ > 0, Δt * τₑₒ / mᵢ / β, 0)

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
                                   v_top_stress,
                                   v_forcing,
                                   model_fields)

    i, j = @index(Global, NTuple)

    uᵢ, vᵢ = velocities
    uₒ, vₒ = ocean_velocities
    h  = thickness
    ℵ  = concentration
    ρᵢ = 917

    vⁿ = rheology.vⁿ

    mᵢ = @inbounds h[i, j, 1] * ℵ[i, j, 1] * ρᵢ
    Δu = ℑxyᶜᶠᵃ(i, j, 1, grid, uₒ) - ℑxyᶜᶠᵃ(i, j, 1, grid, uᵢ)
    Δv = @inbounds vₒ[i, j, 1] - vᵢ[i, j, 1]
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    β = beta_coefficient(rheology, Δt)

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * vₒ - τₑₒ * vᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) / mᵢ
    τₑₒ = 1e-3 * 1020 * Δ𝒰

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τva = v_top_stress[i, j, 1]
    
    @inbounds Gⱽ = ( - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
                     + τva
                     + τₑₒ * vₒ[i, j, 1] # Explicit component of the ice-ocean stress
                     + y_internal_stress_divergence(i, j, grid, rheology) ) / mᵢ

    # make sure we do not have NaNs!
    Gⱽ = ifelse(mᵢ > 0, Gⱽ, 0) 

    # Explicit step
    @inbounds vᵢ[i, j, 1] = (vᵢ[i, j, 1] * (β - 1) + Δt * Gⱽ + vⁿ[i, j, 1]) / β

    # Implicit component of the ice-ocean stress
    τᵢ = ifelse(mᵢ > 0, Δt * τₑₒ / mᵢ / β, 0)

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
