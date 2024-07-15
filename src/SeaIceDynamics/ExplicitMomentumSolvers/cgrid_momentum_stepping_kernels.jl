using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U
using ClimaSeaIce.SeaIceDynamics: Vᵢ

# The ice-ocean stress is treated semi-implicitly 
# i.e:
#
#          Cᴰρₒ
# τₒ =    ------ || uₒ - uⁿ ||   * (uₒ - uⁿ⁺¹)
#           mᵢ
#      |-----------------------|
# τₒ =  τₑ₀ (explicit component)  * Δu   
#

""" stepping the ice u-velocity using a forward leap-frog scheme """
@kernel function _u_cgrid_velocity_step!(velocities, grid, Δt, 
                                         clock,
                                         ocean_velocities,
                                         coriolis,
                                         rheology,
                                         auxiliary_fields,
                                         substeps,
                                         substepping_coefficient,
                                         thickness,
                                         concentration,
                                         ice_density,
                                         ocean_density,
                                         ocean_ice_drag_coefficient,
                                         u_top_stress,
                                         u_forcing,
                                         model_fields)

    i, j = @index(Global, NTuple)

    uᵢ, vᵢ = velocities
    uₒ, vₒ = ocean_velocities
    h  = thickness
    ℵ  = concentration
    ρᵢ = ice_density
    ρₒ = ocean_density
    Cᴰ = ocean_ice_drag_coefficient

    hf = ℑxᴮᶠᶜᶜ(i, j, 1, grid, h) # thickness
    ℵf = ℑxᴮᶠᶜᶜ(i, j, 1, grid, ℵ) # concentration

    # Ice mass (per unit area) interpolated on u points
    mᵢ = hf * ℵf * ρᵢ

    # relative ocean - ice velocities
    Δu = @inbounds uₒ[i, j, 1] - uᵢ[i, j, 1]
    Δv = ℑxyᴮᶠᶜᶜ(i, j, 1, grid, vₒ) - ℑxyᴮᶠᶜᶜ(i, j, 1, grid, vᵢ)

    # relative ocean - ice speed
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    β = get_stepping_coefficients(i, j, 1, grid, substeps, substepping_coefficient)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τuₐ = u_top_stress[i, j, 1] / mᵢ

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * uₒ - τₑₒ * uᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) / mᵢ
    τₑₒ = Cᴰ * ρₒ * Δ𝒰 / mᵢ

    @inbounds Gᵁ = ( - x_f_cross_U(i, j, 1, grid, coriolis, velocities) 
                     + τuₐ
                     + τₑₒ * uₒ[i, j, 1] # Explicit component of the ice-ocean stress
                     + x_internal_stress_divergenceᶠᶜᶜ(i, j, 1, grid, auxiliary_fields, rheology) / mᵢ)

    # make sure we do not have NaNs!                 
    Gᵁ = ifelse(mᵢ > 0, Gᵁ, zero(0)) 
    Gᴿ = rheology_specific_numerical_terms_xᶠᶜᶜ(i, j, 1, grid, rheology, auxiliary_fields, uᵢ)
    
    # Explicit step
    @inbounds uᵢ[i, j, 1] += (Δt * Gᵁ + Gᴿ) / β
    
    # Implicit component of the ice-ocean stress
    τᵢ = ifelse(mᵢ > 0, Δt * τₑₒ / β, zero(grid))

    # Implicit step
    @inbounds uᵢ[i, j, 1] /= (1 + τᵢ) 
end

""" stepping the ice v-velocity using a forward leap-frog scheme """
@kernel function _v_cgrid_velocity_step!(velocities, grid, Δt, 
                                         clock,
                                         ocean_velocities, 
                                         coriolis,
                                         rheology,
                                         auxiliary_fields,
                                         substeps,
                                         substepping_coefficient,
                                         thickness,
                                         concentration,
                                         ice_density,
                                         ocean_density,
                                         ocean_ice_drag_coefficient,
                                         v_top_stress,
                                         v_forcing,
                                         model_fields)

    i, j = @index(Global, NTuple)

    uᵢ, vᵢ = velocities
    uₒ, vₒ = ocean_velocities
    h  = thickness
    ℵ  = concentration
    ρᵢ = ice_density
    ρₒ = ocean_density
    Cᴰ = ocean_ice_drag_coefficient

    hf = ℑyᴮᶜᶠᶜ(i, j, 1, grid, h)
    ℵf = ℑyᴮᶜᶠᶜ(i, j, 1, grid, ℵ)

    # Ice mass interpolated on v points
    mᵢ = hf * ℵf * ρᵢ
    
    # relative ocean - ice velocities
    Δu = ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uₒ) - ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uᵢ)
    Δv = @inbounds vₒ[i, j, 1] - vᵢ[i, j, 1]

    # relative ocean - ice speed
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    β = get_stepping_coefficients(i, j, 1, grid, substeps, substepping_coefficient)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τva = v_top_stress[i, j, 1] / mᵢ 

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * vₒ - τₑₒ * vᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) / mᵢ
    τₑₒ = Cᴰ * ρₒ * Δ𝒰 / mᵢ

    @inbounds Gⱽ = ( - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
                     + τva
                     + τₑₒ * vₒ[i, j, 1] # Explicit component of the ice-ocean stress
                     + y_internal_stress_divergenceᶜᶠᶜ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢ) 

    # make sure we do not have NaNs!
    Gⱽ = ifelse(mᵢ > 0, Gⱽ, zero(0)) 
    Gᴿ = rheology_specific_numerical_terms_yᶜᶠᶜ(i, j, 1, grid, rheology, auxiliary_fields, vᵢ)

    # Explicit step
    @inbounds vᵢ[i, j, 1] += (Δt * Gⱽ + Gᴿ) / β

    # Implicit component of the ice-ocean stress
    τᵢ = ifelse(mᵢ > 0, Δt * τₑₒ / β, zero(0)) 

    # Implicit step
    @inbounds vᵢ[i, j, 1] /= (1 + τᵢ) 
end