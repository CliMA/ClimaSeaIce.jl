using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U, fᶠᶠᵃ

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
                                         immersed_bc,
                                         clock,
                                         ocean_velocities,
                                         coriolis,
                                         rheology,
                                         auxiliary_fields,
                                         substeps,
                                         substepping_coefficient,
                                         ice_thickness,
                                         ice_concentration,
                                         ice_density,
                                         ocean_density,
                                         ocean_ice_drag_coefficient,
                                         u_top_stress,
                                         u_forcing,
                                         model_fields)

    i, j = @index(Global, NTuple)

    uᵢ, vᵢ = velocities.u, velocities.v
    uₒ, vₒ = ocean_velocities.u, ocean_velocities.v
    ûᵢ, v̂ᵢ = auxiliary_fields.û, auxiliary_fields.v̂
    h  = ice_thickness
    ℵ  = ice_concentration
    ρᵢ = ice_density
    ρₒ = ocean_density
    Cᴰ = ocean_ice_drag_coefficient

    # Ice mass (per unit area) interpolated on u points
    mᵢᶠᶜᶜ = ℑxᴮᶠᶜᶜ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)
    mᵢᶜᶠᶜ = ℑyᴮᶜᶠᶜ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

    # relative ocean - ice velocities
    Δuᶠᶜᶜ = @inbounds uₒ[i, j, 1] - uᵢ[i, j, 1]
    Δvᶠᶜᶜ = @inbounds ℑxyᴮᶠᶜᶜ(i, j, 1, grid, vₒ) - v̂ᵢ[i, j, 1]

    Δuᶜᶠᶜ = @inbounds ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uₒ) - ûᵢ[i, j, 1]
    Δvᶜᶠᶜ = @inbounds vₒ[i, j, 1] - vᵢ[i, j, 1]

    # relative ocean - ice speed
    Δ𝒰ᶠᶜᶜ = sqrt(Δuᶠᶜᶜ^2 + Δvᶠᶜᶜ^2)
    Δ𝒰ᶜᶠᶜ = sqrt(Δuᶜᶠᶜ^2 + Δvᶜᶠᶜ^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    βᶠᶜᶜ = ℑxᴮᶠᶜᶜ(i, j, 1, grid, get_stepping_coefficients, substeps, substepping_coefficient)
    βᶜᶠᶜ = ℑyᴮᶜᶠᶜ(i, j, 1, grid, get_stepping_coefficients, substeps, substepping_coefficient)

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
                        + ∂ⱼ_σ₁ⱼᶠᶜᶜ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢᶠᶜᶜ)

    @inbounds Gᵁᶜᶠᶜ = ( + fᶜᶠᶜ(i, j, 1, grid, coriolis) * vᵢ[i, j, 1] 
                        + τuₐᶜᶠᶜ
                        + τₑₒᶜᶠᶜ * ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uₒ) # Explicit component of the ice-ocean stress
                        + ∂ⱼ_σ₁ⱼᶜᶠᶜ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢᶜᶠᶜ)

    # make sure we do not have NaNs!                 
    Gᵁᶠᶜᶜ = ifelse(mᵢᶠᶜᶜ > 0, Gᵁᶠᶜᶜ, zero(grid)) 
    Gᵁᶜᶠᶜ = ifelse(mᵢᶜᶠᶜ > 0, Gᵁᶜᶠᶜ, zero(grid)) 

    Gᴿᶠᶜᶜ = rheology_specific_forcing_xᶠᶜᶜ(i, j, 1, grid, rheology, auxiliary_fields, uᵢ)
    Gᴿᶜᶠᶜ = rheology_specific_forcing_xᶜᶠᶜ(i, j, 1, grid, rheology, auxiliary_fields, ûᵢ)
    
    # Explicit step
    @inbounds uᵢ[i, j, 1] += (Δt * Gᵁᶠᶜᶜ + Gᴿᶠᶜᶜ) / βᶠᶜᶜ
    @inbounds ûᵢ[i, j, 1] += (Δt * Gᵁᶜᶠᶜ + Gᴿᶜᶠᶜ) / βᶜᶠᶜ
    
    # Implicit component of the ice-ocean stress
    τᵢᶠᶜᶜ = ifelse(mᵢᶠᶜᶜ > 0, Δt * τₑₒᶠᶜᶜ / βᶠᶜᶜ, zero(grid))
    τᵢᶜᶠᶜ = ifelse(mᵢᶜᶠᶜ > 0, Δt * τₑₒᶜᶠᶜ / βᶜᶠᶜ, zero(grid))

    # Implicit step
    @inbounds uᵢ[i, j, 1] /= (1 + τᵢᶠᶜᶜ) 
    @inbounds ûᵢ[i, j, 1] /= (1 + τᵢᶜᶠᶜ) 
end

""" stepping the ice v-velocity using a forward leap-frog scheme """
@kernel function _v_egrid_velocity_step!(velocities, grid, Δt, 
                                         immersed_bc,
                                         clock,
                                         ocean_velocities, 
                                         coriolis,
                                         rheology,
                                         auxiliary_fields,
                                         substeps,
                                         substepping_coefficient,
                                         ice_thickness,
                                         ice_concentration,
                                         ice_density,
                                         ocean_density,
                                         ocean_ice_drag_coefficient,
                                         v_top_stress,
                                         v_forcing,
                                         model_fields)

    i, j = @index(Global, NTuple)

    uᵢ, vᵢ = velocities.u, velocities.v
    uₒ, vₒ = ocean_velocities.u, ocean_velocities.v
    ûᵢ, v̂ᵢ = auxiliary_fields.û, auxiliary_fields.v̂
    h  = ice_thickness
    ℵ  = ice_concentration
    ρᵢ = ice_density
    ρₒ = ocean_density
    Cᴰ = ocean_ice_drag_coefficient

    # Ice mass (per unit area) interpolated on u points
    mᵢᶠᶜᶜ = ℑxᴮᶠᶜᶜ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)
    mᵢᶜᶠᶜ = ℑyᴮᶜᶠᶜ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)
    
    # relative ocean - ice velocities
    Δuᶠᶜᶜ = @inbounds uₒ[i, j, 1] - uᵢ[i, j, 1]
    Δvᶠᶜᶜ = @inbounds ℑxyᴮᶠᶜᶜ(i, j, 1, grid, vₒ) - v̂ᵢ[i, j, 1]

    Δuᶜᶠᶜ = @inbounds ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uₒ) - ûᵢ[i, j, 1]
    Δvᶜᶠᶜ = @inbounds vₒ[i, j, 1] - vᵢ[i, j, 1]

    # relative ocean - ice speed
    Δ𝒰ᶠᶜᶜ = sqrt(Δuᶠᶜᶜ^2 + Δvᶠᶜᶜ^2)
    Δ𝒰ᶜᶠᶜ = sqrt(Δuᶜᶠᶜ^2 + Δvᶜᶠᶜ^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    βᶠᶜᶜ = ℑxᴮᶠᶜᶜ(i, j, 1, grid, get_stepping_coefficients, substeps, substepping_coefficient)
    βᶜᶠᶜ = ℑyᴮᶜᶠᶜ(i, j, 1, grid, get_stepping_coefficients, substeps, substepping_coefficient)

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
                        + ∂ⱼ_σ₂ⱼᶜᶠᶜ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢᶜᶠᶜ) 

    @inbounds Gⱽᶠᶜᶜ = ( - fᶠᶜᶜ(i, j, 1, grid, coriolis) * uᵢ[i, j, 1] 
                        + τvₐᶠᶜᶜ
                        + τₑₒᶠᶜᶜ * vₒ[i, j, 1] # Explicit component of the ice-ocean stress
                        + ∂ⱼ_σ₂ⱼᶠᶜᶜ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢᶠᶜᶜ) 

    # make sure we do not have NaNs!
    Gⱽᶜᶠᶜ = ifelse(mᵢᶜᶠᶜ > 0, Gⱽᶜᶠᶜ, zero(grid)) 
    Gⱽᶠᶜᶜ = ifelse(mᵢᶠᶜᶜ > 0, Gⱽᶠᶜᶜ, zero(grid)) 

    Gᴿᶜᶠᶜ = rheology_specific_forcing_yᶜᶠᶜ(i, j, 1, grid, rheology, auxiliary_fields, vᵢ)
    Gᴿᶠᶜᶜ = rheology_specific_forcing_yᶠᶜᶜ(i, j, 1, grid, rheology, auxiliary_fields, v̂ᵢ)

    # Explicit step
    @inbounds vᵢ[i, j, 1] += (Δt * Gⱽᶜᶠᶜ + Gᴿᶜᶠᶜ) / βᶜᶠᶜ
    @inbounds v̂ᵢ[i, j, 1] += (Δt * Gⱽᶠᶜᶜ + Gᴿᶠᶜᶜ) / βᶠᶜᶜ

    # Implicit component of the ice-ocean stress
    τᵢᶜᶠᶜ = ifelse(mᵢᶜᶠᶜ > 0, Δt * τₑₒᶜᶠᶜ / βᶜᶠᶜ, zero(grid)) 
    τᵢᶠᶜᶜ = ifelse(mᵢᶠᶜᶜ > 0, Δt * τₑₒᶠᶜᶜ / βᶠᶜᶜ, zero(grid)) 

    # Implicit step
    @inbounds vᵢ[i, j, 1] /= (1 + τᵢᶜᶠᶜ) 
    @inbounds v̂ᵢ[i, j, 1] /= (1 + τᵢᶠᶜᶜ) 
end