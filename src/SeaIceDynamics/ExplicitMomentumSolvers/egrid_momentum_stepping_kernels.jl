using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U, fᶠᶠᵃ
using ClimaSeaIce.SeaIceDynamics: Vᵢ


@inline fᶠᶜᶜ(i, j, k, grid, coriolis) = ℑyᴮᶠᶜᶜ(i, j, k, grid, fᶠᶠᵃ, coriolis)
@inline fᶜᶠᶜ(i, j, k, grid, coriolis) = ℑxᴮᶜᶠᶜ(i, j, k, grid, fᶠᶠᵃ, coriolis)

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
@kernel function _u_egrid_velocity_step!(velocities, grid, Δt, 
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
    ûᵢ, v̂ᵢ = auxiliary_fields.û, auxiliary_fields.v̂
    h  = thickness
    ℵ  = concentration
    ρᵢ = ice_density
    ρₒ = ocean_density
    Cᴰ = ocean_ice_drag_coefficient

    hᶠᶜᶜ = ℑxᴮᶠᶜᶜ(i, j, 1, grid, h) # thickness
    ℵᶠᶜᶜ = ℑxᴮᶠᶜᶜ(i, j, 1, grid, ℵ) # concentration

    hᶜᶠᶜ = ℑyᴮᶜᶠᶜ(i, j, 1, grid, h) # thickness
    ℵᶜᶠᶜ = ℑyᴮᶜᶠᶜ(i, j, 1, grid, ℵ) # concentration

    # Ice mass (per unit area) interpolated on u points
    mᵢᶠᶜᶜ = hᶠᶜᶜ * ℵᶠᶜᶜ * ρᵢ
    mᵢᶜᶠᶜ = hᶜᶠᶜ * ℵᶜᶠᶜ * ρᵢ

    # relative ocean - ice velocities
    Δuᶠᶜᶜ = @inbounds uₒ[i, j, 1] - uᵢ[i, j, 1]
    Δvᶠᶜᶜ = @inbounds ℑxyᴮᶠᶜᶜ(i, j, 1, grid, vₒ) - v̂ᵢ[i, j, 1]

    Δuᶜᶠᶜ = @inbounds ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uₒ) - ûᵢ[i, j, 1]
    Δvᶜᶠᶜ = @inbounds vₒ[i, j, 1] - vᵢ[i, j, 1]

    # relative ocean - ice speed
    Δ𝒰ᶠᶜᶜ = sqrt(Δuᶠᶜᶜ^2 + Δvᶠᶜᶜ^2)
    Δ𝒰ᶜᶠᶜ = sqrt(Δuᶜᶠᶜ^2 + Δvᶜᶠᶜ^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    β = get_stepping_coefficients(i, j, 1, grid, substeps, substepping_coefficient)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τuₐᶠᶜᶜ = u_top_stress[i, j, 1] / mᵢᶠᶜᶜ
    @inbounds τuₐᶜᶠᶜ = ℑxyᴮᶜᶠᶜ(i, j, 1, grid, u_top_stress) / mᵢᶜᶠᶜ

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * uₒ - τₑₒ * uᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) / mᵢ
    τₑₒᶠᶜᶜ = Cᴰ * ρₒ * Δ𝒰ᶠᶜᶜ / mᵢᶠᶜᶜ
    τₑₒᶜᶠᶜ = Cᴰ * ρₒ * Δ𝒰ᶜᶠᶜ / mᵢᶜᶠᶜ

    @inbounds Gᵁᶠᶜᶜ = ( + fᶠᶜᶜ(i, j, 1, grid, coriolis) * v̂ᵢ[i, j, 1] 
                        + τuₐᶠᶜᶜ
                        + τₑₒᶠᶜᶜ * uₒ[i, j, 1] # Explicit component of the ice-ocean stress
                        + x_internal_stress_divergenceᶠᶜᶜ(i, j, 1, grid, auxiliary_fields, rheology) / mᵢᶠᶜᶜ)

    @inbounds Gᵁᶜᶠᶜ = ( + fᶜᶠᶜ(i, j, 1, grid, coriolis) * vᵢ[i, j, 1] 
                        + τuₐᶜᶠᶜ
                        + τₑₒᶜᶠᶜ * ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uₒ) # Explicit component of the ice-ocean stress
                        + x_internal_stress_divergenceᶜᶠᶜ(i, j, 1, grid, auxiliary_fields, rheology) / mᵢᶜᶠᶜ)

    # make sure we do not have NaNs!                 
    Gᵁᶠᶜᶜ = ifelse(mᵢᶠᶜᶜ > 0, Gᵁᶠᶜᶜ, zero(0)) 
    Gᵁᶜᶠᶜ = ifelse(mᵢᶜᶠᶜ > 0, Gᵁᶜᶠᶜ, zero(0)) 
    Gᴿᶠᶜᶜ = rheology_specific_numerical_terms_xᶠᶜᶜ(i, j, 1, grid, rheology, auxiliary_fields, uᵢ)
    Gᴿᶜᶠᶜ = rheology_specific_numerical_terms_xᶜᶠᶜ(i, j, 1, grid, rheology, auxiliary_fields, ûᵢ)
    
    # Explicit step
    @inbounds uᵢ[i, j, 1] += (Δt * Gᵁᶠᶜᶜ + Gᴿᶠᶜᶜ) / β
    @inbounds ûᵢ[i, j, 1] += (Δt * Gᵁᶜᶠᶜ + Gᴿᶜᶠᶜ) / β
    
    # Implicit component of the ice-ocean stress
    τᵢᶠᶜᶜ = ifelse(mᵢ > 0, Δt * τₑₒᶠᶜᶜ / β, zero(grid))
    τᵢᶜᶠᶜ = ifelse(mᵢ > 0, Δt * τₑₒᶜᶠᶜ / β, zero(grid))

    # Implicit step
    @inbounds uᵢ[i, j, 1] /= (1 + τᵢᶠᶜᶜ) 
    @inbounds ûᵢ[i, j, 1] /= (1 + τᵢᶜᶠᶜ) 
end

""" stepping the ice v-velocity using a forward leap-frog scheme """
@kernel function _v_egrid_velocity_step!(velocities, grid, Δt, 
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

    hᶠᶜᶜ = ℑxᴮᶠᶜᶜ(i, j, 1, grid, h) # thickness
    ℵᶠᶜᶜ = ℑxᴮᶠᶜᶜ(i, j, 1, grid, ℵ) # concentration

    hᶜᶠᶜ = ℑyᴮᶜᶠᶜ(i, j, 1, grid, h) # thickness
    ℵᶜᶠᶜ = ℑyᴮᶜᶠᶜ(i, j, 1, grid, ℵ) # concentration

    # Ice mass (per unit area) interpolated on u points
    mᵢᶠᶜᶜ = hᶠᶜᶜ * ℵᶠᶜᶜ * ρᵢ
    mᵢᶜᶠᶜ = hᶜᶠᶜ * ℵᶜᶠᶜ * ρᵢ

    # relative ocean - ice velocities
    Δuᶠᶜᶜ = @inbounds uₒ[i, j, 1] - uᵢ[i, j, 1]
    Δvᶠᶜᶜ = @inbounds ℑxyᴮᶠᶜᶜ(i, j, 1, grid, vₒ) - v̂ᵢ[i, j, 1]

    Δuᶜᶠᶜ = @inbounds ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uₒ) - ûᵢ[i, j, 1]
    Δvᶜᶠᶜ = @inbounds vₒ[i, j, 1] - vᵢ[i, j, 1]

    # relative ocean - ice speed
    Δ𝒰ᶠᶜᶜ = sqrt(Δuᶠᶜᶜ^2 + Δvᶠᶜᶜ^2)
    Δ𝒰ᶜᶠᶜ = sqrt(Δuᶜᶠᶜ^2 + Δvᶜᶠᶜ^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    β = get_stepping_coefficients(i, j, 1, grid, substeps, substepping_coefficient)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τvₐᶠᶜᶜ = v_top_stress[i, j, 1] / mᵢᶠᶜᶜ
    @inbounds τvₐᶜᶠᶜ = ℑxyᴮᶜᶠᶜ(i, j, 1, grid, v_top_stress) / mᵢᶜᶠᶜ

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * uₒ - τₑₒ * uᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) / mᵢ
    τₑₒᶠᶜᶜ = Cᴰ * ρₒ * Δ𝒰ᶠᶜᶜ / mᵢᶠᶜᶜ
    τₑₒᶜᶠᶜ = Cᴰ * ρₒ * Δ𝒰ᶜᶠᶜ / mᵢᶜᶠᶜ

    @inbounds Gⱽᶜᶠᶜ = ( - fᶜᶠᶜ(i, j, 1, grid, coriolis) * ûᵢ[i, j, 1] 
                        + τvₐᶜᶠᶜ
                        + τₑₒᶜᶠᶜ * vₒ[i, j, 1] # Explicit component of the ice-ocean stress
                        + y_internal_stress_divergenceᶜᶠᶜ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢᶜᶠᶜ) 

    @inbounds Gⱽᶠᶜᶜ = ( - fᶠᶜᶜ(i, j, 1, grid, coriolis) * ûᵢ[i, j, 1] 
                        + τvₐᶠᶜᶜ
                        + τₑₒᶠᶜᶜ * vₒ[i, j, 1] # Explicit component of the ice-ocean stress
                        + y_internal_stress_divergenceᶠᶜᶜ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢᶠᶜᶜ) 

    # make sure we do not have NaNs!
    Gⱽᶜᶠᶜ = ifelse(mᵢᶜᶠᶜ > 0, Gⱽᶜᶠᶜ, zero(0)) 
    Gⱽᶠᶜᶜ = ifelse(mᵢᶠᶜᶜ > 0, Gⱽᶠᶜᶜ, zero(0)) 
    Gᴿᶜᶠᶜ = rheology_specific_numerical_terms_yᶜᶠᶜ(i, j, 1, grid, rheology, auxiliary_fields, vᵢ)
    Gᴿᶠᶜᶜ = rheology_specific_numerical_terms_yᶠᶜᶜ(i, j, 1, grid, rheology, auxiliary_fields, v̂ᵢ)

    # Explicit step
    @inbounds vᵢ[i, j, 1] += (Δt * Gⱽᶜᶠᶜ + Gᴿᶜᶠᶜ) / β
    @inbounds v̂ᵢ[i, j, 1] += (Δt * Gⱽᶠᶜᶜ + Gᴿᶠᶜᶜ) / β

    # Implicit component of the ice-ocean stress
    τᵢᶜᶠᶜ = ifelse(mᵢ > 0, Δt * τₑₒᶜᶠᶜ / β, zero(0)) 
    τᵢᶠᶜᶜ = ifelse(mᵢ > 0, Δt * τₑₒᶠᶜᶜ / β, zero(0)) 

    # Implicit step
    @inbounds vᵢ[i, j, 1] /= (1 + τᵢᶜᶠᶜ) 
    @inbounds v̂ᵢ[i, j, 1] /= (1 + τᵢᶠᶜᶜ) 
end