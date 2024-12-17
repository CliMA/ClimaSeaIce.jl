using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U
using Oceananigans.ImmersedBoundaries: active_linear_index_to_tuple

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
@kernel function _u_velocity_step!(velocities, grid, ::Nothing, args...)
    i, j = @index(Global, NTuple)
    u_velocity_step!(i, j, velocities, grid, args...)
end

@kernel function _v_velocity_step!(velocities, grid, ::Nothing, args...)
    i, j = @index(Global, NTuple)
    v_velocity_step!(i, j, velocities, grid, args...)
end

@kernel function _u_velocity_step!(velocities, grid, active_surface_map, args...)
    idx = @index(Global, Linear)
    i, j = active_linear_index_to_tuple(idx, active_surface_map)
    u_velocity_step!(i, j, velocities, grid, args...)
end

@kernel function _v_velocity_step!(velocities, grid, active_surface_map, args...)
    idx = @index(Global, Linear)
    i, j = active_linear_index_to_tuple(idx, active_surface_map)
    v_velocity_step!(i, j, velocities, grid, args...)
end

""" stepping the ice u-velocity using a forward leap-frog scheme """
@inline function u_velocity_step!(i, j, 
                                  velocities, grid, Δt, 
                                  immersed_bc,
                                  clock,
                                  ocean_velocities,
                                  ocean_free_surface,
                                  coriolis,
                                  rheology,
                                  auxiliary_fields,
                                  substeps,
                                  ice_thickness,
                                  ice_concentration,
                                  ice_density,
                                  ocean_ice_drag_coefficient,
                                  gravitational_acceleration,
                                  u_top_stress,
                                  u_forcing)

    uᵢ = velocities.u
    vᵢ = velocities.v
    uₒ = ocean_velocities.u
    vₒ = ocean_velocities.v
    ηₒ = ocean_free_surface
    h  = ice_thickness
    ℵ  = ice_concentration
    ρᵢ = ice_density
    Cᴰ = ocean_ice_drag_coefficient
    g  = gravitational_acceleration

    fields = merge(auxiliary_fields, velocities, (h, ℵ))

    # Ice mass (per unit area) interpolated on u points
    mᵢ = ℑxᶠᵃᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

    # relative ocean - ice velocities
    Δu = @inbounds uₒ[i, j, 1] - uᵢ[i, j, 1]
    Δv = ℑxyᶠᶜᵃ(i, j, 1, grid, vₒ) 
       - ℑxyᶠᶜᵃ(i, j, 1, grid, vᵢ)

    # relative ocean - ice speed
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    β = ℑxᶠᵃᵃ(i, j, 1, grid, rheology_substeps, rheology, substeps, auxiliary_fields)

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
                     + g * ∂xᶠᶜᶜ(i, j, 1, grid, ηₒ)
                     + ∂ⱼ_σ₁ⱼ(i, j, 1, grid, rheology, fields) / mᵢ)

    # make sure we do not have NaNs!                 
    Gᵁ = ifelse(mᵢ > 0, Gᵁ, zero(grid)) 
    Gᴿ = rheology_specific_forcing_x(i, j, 1, grid, rheology, fields)
    
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
                                  ocean_free_surface,
                                  coriolis,
                                  rheology,
                                  auxiliary_fields,
                                  substeps,
                                  ice_thickness,
                                  ice_concentration,
                                  ice_density,
                                  ocean_ice_drag_coefficient,
                                  gravitational_acceleration,
                                  v_top_stress,
                                  v_forcing)

    uᵢ = velocities.u
    vᵢ = velocities.v
    uₒ = ocean_velocities.u
    vₒ = ocean_velocities.v
    ηₒ = ocean_free_surface
    h  = ice_thickness
    ℵ  = ice_concentration
    ρᵢ = ice_density
    Cᴰ = ocean_ice_drag_coefficient
    g  = gravitational_acceleration

    fields = merge(auxiliary_fields, velocities, (h, ℵ))

    # Ice mass (per unit area) interpolated on u points
    mᵢ = ℑyᵃᶠᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

    # relative ocean - ice velocities
    Δu = ℑxyᶜᶠᵃ(i, j, 1, grid, uₒ) 
       - ℑxyᶜᶠᵃ(i, j, 1, grid, uᵢ)

    Δv = @inbounds vₒ[i, j, 1] - vᵢ[i, j, 1]

    # relative ocean - ice speed
    Δ𝒰 = sqrt(Δu^2 + Δv^2)
    
    # Coefficient for substepping momentum (depends on the particular substepping formulation)
    β = ℑyᵃᶠᵃ(i, j, 1, grid, rheology_substeps, rheology, substeps, fields)

    # The atmosphere - ice stress is prescribed at each time step
    # (i.e. it only depends on wind speed)
    @inbounds τvₐ = v_top_stress[i, j, 1] / mᵢ 

    # The ocean - ice stress is computed semi-implicitly as
    # τₒ = τₑₒ * vₒ - τₑₒ * vᵢⁿ⁺¹ 
    # where τₑₒ = (Cᴰ ρₒ Δ𝒰ⁿ) / mᵢ
    τₑₒ = Cᴰ * Δ𝒰 / mᵢ

    @inbounds Gⱽ =  - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
                    + τvₐ
                    + g * ∂yᶜᶠᶜ(i, j, 1, grid, ηₒ)
                    + τₑₒ * vₒ[i, j, 1] # Explicit component of the ice-ocean stress
                    + ∂ⱼ_σ₂ⱼ(i, j, 1, grid, rheology, auxiliary_fields) / mᵢ

    # make sure we do not have NaNs!
    Gⱽ = ifelse(mᵢ > 0, Gⱽ, zero(grid)) 
    Gᴿ = rheology_specific_forcing_y(i, j, 1, grid, rheology, fields)

    # Explicit step
    @inbounds vᵢ[i, j, 1] += (Δt * Gⱽ + Gᴿ) / β

    # Implicit component of the ice-ocean stress
    τᵢ = ifelse(mᵢ > 0, Δt * τₑₒ / β, zero(grid)) 

    # Implicit step
    @inbounds vᵢ[i, j, 1] /= (1 + τᵢ) 
end