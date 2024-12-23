using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U
using Oceananigans.ImmersedBoundaries: active_linear_index_to_tuple

"""compute ice u-velocity tendencies"""
@inline function u_velocity_tendency(i, j, grid,
                                     clock,
                                     velocities,
                                     coriolis,
                                     rheology,
                                     auxiliary_fields,
                                     ice_thickness,
                                     ice_concentration,
                                     ice_density,
                                     u_top_stress,
                                     u_bottom_stress)

   h  = ice_thickness
   ℵ  = ice_concentration
   ρᵢ = ice_density

   fields = merge(auxiliary_fields, velocities, (; h, ℵ))

   # Ice mass (per unit area) interpolated on u points
   mᵢ = ℑxᶠᵃᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

   @inbounds Gᵁ = ( - x_f_cross_U(i, j, 1, grid, coriolis, velocities) 
                    + τx(i, j, 1, grid, u_top_stress, velocities, mᵢ)
                    + τx(i, j, 1, grid, u_bottom_stress, velocities, mᵢ)
                    + ∂ⱼ_σ₁ⱼ(i, j, 1, grid, rheology, clock, fields) / mᵢ)

   return ifelse(mᵢ ≤ 0, zero(grid), Gᵁ) 
end

"""compute ice v-velocity tendencies"""
@inline function v_velocity_tendency!(i, j, grid,
                                      clock,
                                      velocities,
                                      coriolis,
                                      rheology,
                                      auxiliary_fields,
                                      ice_thickness,
                                      ice_concentration,
                                      ice_density,
                                      v_top_stress,
                                      v_bottom_stress)

   h  = ice_thickness
   ℵ  = ice_concentration
   ρᵢ = ice_density

   fields = merge(auxiliary_fields, velocities, (; h, ℵ))

   # Ice mass (per unit area) interpolated on u points
   mᵢ = ℑyᵃᶠᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

   @inbounds Gⱽ = ( - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
                    + τy(i, j, 1, grid, v_top_stress, velocities, mᵢ)
                    + τy(i, j, 1, grid, v_bottom_stress, velocities, mᵢ)
                    + ∂ⱼ_σ₂ⱼ(i, j, 1, grid, rheology, clock, fields) / mᵢ)

   return ifelse(mᵢ ≤ 0, zero(grid), Gⱽ) 
end

@inline τx(i, j, k, grid, stress::Nothing, velocities, mᵢ) = zero(grid)
@inline τy(i, j, k, grid, stress::Nothing, velocities, mᵢ) = zero(grid)

@inline τx(i, j, k, grid, stress::Number, velocities, mᵢ) = stress / mᵢ
@inline τy(i, j, k, grid, stress::Number, velocities, mᵢ) = stress / mᵢ

@inline τx(i, j, k, grid, stress::AbstractArray, velocities, mᵢ) =  @inbounds stress[i, j, k] / mᵢ
@inline τy(i, j, k, grid, stress::AbstractArray, velocities, mᵢ) =  @inbounds stress[i, j, k] / mᵢ
