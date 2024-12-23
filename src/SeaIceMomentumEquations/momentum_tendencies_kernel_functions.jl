using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U
using Oceananigans.ImmersedBoundaries: active_linear_index_to_tuple

""" stepping the ice u-velocity using a forward leap-frog scheme """
@inline function u_velocity_tendency(i, j, grid,
                                     clock,
                                     velocities,
                                     ocean_free_surface,
                                     coriolis,
                                     rheology,
                                     auxiliary_fields,
                                     ice_thickness,
                                     ice_concentration,
                                     ice_density,
                                     gravitational_acceleration,
                                     u_top_stress,
                                     u_bottom_stress)

   ηₒ = ocean_free_surface
   h  = ice_thickness
   ℵ  = ice_concentration
   ρᵢ = ice_density
   g  = gravitational_acceleration

   fields = merge(auxiliary_fields, velocities, (; h, ℵ))

   # Ice mass (per unit area) interpolated on u points
   mᵢ = ℑxᶠᵃᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

   @inbounds Gᵁ = ( - x_f_cross_U(i, j, 1, grid, coriolis, velocities) 
                    + τ_atmosphere_x(i, j, 1, grid, u_top_stress, velocities, mᵢ)
                    + τ_ocean_x(i, j, 1, grid, u_bottom_stress, velocities, mᵢ)
                    + g * ∂xᶠᶜᶜ(i, j, 1, grid, ηₒ)
                    + ∂ⱼ_σ₁ⱼ(i, j, 1, grid, rheology, clock, fields) / mᵢ)

   return ifelse(mᵢ ≤ 0, zero(grid), Gᵁ) 
end

""" stepping the ice v-velocity using a forward leap-frog scheme """
@inline function v_velocity_tendency!(i, j, grid,
                                      clock,
                                      velocities,
                                      ocean_free_surface,
                                      coriolis,
                                      rheology,
                                      auxiliary_fields,
                                      ice_thickness,
                                      ice_concentration,
                                      ice_density,
                                      gravitational_acceleration,
                                      v_top_stress,
                                      v_bottom_stress)

   ηₒ = ocean_free_surface
   h  = ice_thickness
   ℵ  = ice_concentration
   ρᵢ = ice_density
   g  = gravitational_acceleration

   fields = merge(auxiliary_fields, velocities, (; h, ℵ))

   # Ice mass (per unit area) interpolated on u points
   mᵢ = ℑyᵃᶠᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

   @inbounds Gⱽ = ( - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
                    + τ_atmosphere_y(i, j, 1, grid, v_top_stress, velocities, mᵢ)
                    + τ_ocean_y(i, j, 1, grid, v_bottom_stress, velocities, mᵢ)
                    + g * ∂yᶜᶠᶜ(i, j, 1, grid, ηₒ)
                    + ∂ⱼ_σ₂ⱼ(i, j, 1, grid, rheology, clock, fields) / mᵢ)

   return ifelse(mᵢ ≤ 0, zero(grid), Gⱽ) 
end