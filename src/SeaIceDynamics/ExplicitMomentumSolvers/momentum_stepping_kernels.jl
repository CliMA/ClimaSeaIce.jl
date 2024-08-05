using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U
using Oceananigans.ImmersedBoundaries: active_linear_index_to_tuple, ActiveZColumnsIBG, ZColumnMap
# The ice-ocean stress is treated semi-implicitly 
# i.e:
#
#          Cᴰρₒ
# τₒ =    ------ || uₒ - uⁿ ||   * (uₒ - uⁿ⁺¹)
#           mᵢ
#      |-----------------------|
# τₒ =  τₑ₀ (explicit component)  * Δu   
#

# Make sure we do not compute inside the immersed boundary
@kernel function _u_velocity_step!(velocities, grid, args...)
    i, j = @index(Global, NTuple)
    u_velocity_step!(i, j, velocities, grid, args...)
end

@kernel function _v_velocity_step!(velocities, grid, args...)
    i, j = @index(Global, NTuple)
    v_velocity_step!(i, j, velocities, grid, args...)
end

@kernel function _u_velocity_step!(velocities, grid::ActiveZColumnsIBG, args...)
    idx = @index(Global, Linear)
    i, j = active_linear_index_to_tuple(idx, ZColumnMap(), grid)
    u_velocity_step!(i, j, velocities, grid, args...)
end

@kernel function _v_velocity_step!(velocities, grid::ActiveZColumnsIBG, args...)
    idx = @index(Global, Linear)
    i, j = active_linear_index_to_tuple(idx, ZColumnMap(), grid)
    v_velocity_step!(i, j, velocities, grid, args...)
end

""" stepping the ice u-velocity using a forward leap-frog scheme """
@inline function u_velocity_step!(i, j, 
                                  velocities, grid, Δt, 
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
                                  ocean_ice_drag_coefficient,
                                  u_top_stress,
                                  u_forcing,
                                  model_fields)

    uᵢ = velocities.u
    vᵢ = velocities.v
    uₒ = ocean_velocities.u
    vₒ = ocean_velocities.v
    h  = ice_thickness
    ℵ  = ice_concentration
    ρᵢ = ice_density
    Cᴰ = ocean_ice_drag_coefficient

    # Ice mass (per unit area) interpolated on u points
    mᵢ = ℑxᴮᶠᶜᶜ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

    # relative ocean - ice velocities
    Δu = @inbounds uₒ[i, j, 1] - uᵢ[i, j, 1]
    Δv = ℑxyᴮᶠᶜᶜ(i, j, 1, grid, vₒ) - ℑxyᴮᶠᶜᶜ(i, j, 1, grid, vᵢ)

    # relative ocean - ice speed
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    β = ℑxᴮᶠᶜᶜ(i, j, 1, grid, get_stepping_coefficients, substeps, substepping_coefficient)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τuₐ = u_top_stress[i, j, 1] / mᵢ

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * uₒ - τₑₒ * uᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) / mᵢ
    τₑₒ = Cᴰ * Δ𝒰 / mᵢ

    @inbounds Gᵁ = ( - x_f_cross_U(i, j, 1, grid, coriolis, velocities) 
                     + τuₐ
                     + τₑₒ * uₒ[i, j, 1] # Explicit component of the ice-ocean stress
                     + ∂ⱼ_σ₁ⱼ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢ)

    # make sure we do not have NaNs!                 
    Gᵁ = ifelse(mᵢ > 0, Gᵁ, zero(grid)) 
    Gᴿ = rheology_specific_forcing_x(i, j, 1, grid, rheology, auxiliary_fields, uᵢ)
    
    # Explicit step
    @inbounds uᵢ[i, j, 1] += (Δt * Gᵁ + Gᴿ) / β
    
    # Implicit component of the ice-ocean stress
    τᵢ = ifelse(mᵢ > 0, Δt * τₑₒ / β, zero(grid))

    # Implicit step
    @inbounds uᵢ[i, j, 1] /= (1 + τᵢ) 
end

""" stepping the ice v-velocity using a forward leap-frog scheme """
@inline function v_velocity_step!(i, j, 
                                  velocities, grid, Δt, 
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
                                  ocean_ice_drag_coefficient,
                                  v_top_stress,
                                  v_forcing,
                                  model_fields)

    uᵢ = velocities.u
    vᵢ = velocities.v
    uₒ = ocean_velocities.u
    vₒ = ocean_velocities.v
    h  = ice_thickness
    ℵ  = ice_concentration
    ρᵢ = ice_density
    Cᴰ = ocean_ice_drag_coefficient

    # Ice mass (per unit area) interpolated on u points
    mᵢ = ℑyᴮᶜᶠᶜ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

    # relative ocean - ice velocities
    Δu = ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uₒ) - ℑxyᴮᶜᶠᶜ(i, j, 1, grid, uᵢ)
    Δv = @inbounds vₒ[i, j, 1] - vᵢ[i, j, 1]

    # relative ocean - ice speed
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    β = ℑyᴮᶜᶠᶜ(i, j, 1, grid, get_stepping_coefficients, substeps, substepping_coefficient)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τvₐ = v_top_stress[i, j, 1] / mᵢ 

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * vₒ - τₑₒ * vᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) / mᵢ
    τₑₒ = Cᴰ * Δ𝒰 / mᵢ

    @inbounds Gⱽ = ( - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
                     + τvₐ
                     + τₑₒ * vₒ[i, j, 1] # Explicit component of the ice-ocean stress
                     + ∂ⱼ_σ₂ⱼ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢ) 

    # make sure we do not have NaNs!
    Gⱽ = ifelse(mᵢ > 0, Gⱽ, zero(grid)) 
    Gᴿ = rheology_specific_forcing_y(i, j, 1, grid, rheology, auxiliary_fields, vᵢ)

    # Explicit step
    @inbounds vᵢ[i, j, 1] += (Δt * Gⱽ + Gᴿ) / β

    # Implicit component of the ice-ocean stress
    τᵢ = ifelse(mᵢ > 0, Δt * τₑₒ / β, zero(grid)) 

    # Implicit step
    @inbounds vᵢ[i, j, 1] /= (1 + τᵢ) 
end