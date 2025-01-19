using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U
using Oceananigans.ImmersedBoundaries: active_linear_index_to_tuple

"""compute ice u-velocity tendencies"""
@inline function u_velocity_tendency(i, j, grid, Δt,
                                     rheology,
                                     auxiliary_fields,
                                     clock,
                                     velocities,
                                     coriolis,
                                     ice_thickness,
                                     ice_concentration,
                                     ice_density,
                                     u_top_stress,
                                     u_bottom_stress)

   h = ice_thickness
   ℵ = ice_concentration
   ρ = ice_density

   fields = merge(auxiliary_fields, velocities, (; h, ℵ))

   # Ice mass (per unit area) interpolated on u points
   mᵢ = ℑxᶠᵃᵃ(i, j, 1, grid, ice_mass, h, ρ)
   ℵᵢ = ℑyᵃᶠᵃ(i, j, 1, grid, ℵ)

   Gᵁ = ( - x_f_cross_U(i, j, 1, grid, coriolis, velocities) 
          + explicit_τx(i, j, 1, grid, u_top_stress, clock, fields) / mᵢ * ℵᵢ
          + explicit_τx(i, j, 1, grid, u_bottom_stress, clock, fields) / mᵢ * ℵᵢ
          + ∂ⱼ_σ₁ⱼ(i, j, 1, grid, rheology, clock, fields, Δt) / mᵢ
          # sum of user defined forcing and possibly other forcing terms that are rheology-dependent 
          + sum_of_forcing_x(i, j, 1, grid, rheology, u_forcing, fields, Δt)) 

   # Implicit part of the stress that depends linearly on the velocity
   τᵢ = ( implicit_τx_coefficient(i, j, 1, grid, u_bottom_stress, clock, fields) / mᵢ * ℵᵢ
        + implicit_τx_coefficient(i, j, 1, grid, u_top_stress, clock, fields) / mᵢ * ℵᵢ )

   Gᵁ = ifelse(mᵢ ≤ 0, zero(grid), Gᵁ)
   τᵢ = ifelse(mᵢ ≤ 0, zero(grid), τᵢ)

   return τᵢ, Gᵁ
end

"""compute ice v-velocity tendencies"""
@inline function v_velocity_tendency(i, j, grid, Δt,
                                     rheology,
                                     auxiliary_fields,
                                     clock,
                                     velocities,
                                     coriolis,
                                     ice_thickness,
                                     ice_concentration,
                                     ice_density,
                                     v_top_stress,
                                     v_bottom_stress,
                                     v_forcing)

   h = ice_thickness
   ℵ = ice_concentration
   ρ = ice_density

   fields = merge(auxiliary_fields, velocities, (; h, ℵ))

   # Ice mass (per unit area) interpolated on v points
   mᵢ = ℑyᵃᶠᵃ(i, j, 1, grid, ice_mass, h, ρ)
   ℵᵢ = ℑyᵃᶠᵃ(i, j, 1, grid, ℵ)

    Gⱽ = ( - y_f_cross_U(i, j, 1, grid, coriolis, velocities)
           + explicit_τy(i, j, 1, grid, v_top_stress, clock, fields) / mᵢ * ℵᵢ
           + explicit_τy(i, j, 1, grid, v_bottom_stress, clock, fields) / mᵢ * ℵᵢ
           + ∂ⱼ_σ₂ⱼ(i, j, 1, grid, rheology, clock, fields, Δt) / mᵢ 
           # sum of user defined forcing and possibly other forcing terms that are rheology-dependent 
          + sum_of_forcing_y(i, j, 1, grid, rheology, v_forcing, fields, Δt))

   # Implicit part of the stress that depends linearly on the velocity
   τᵢ = ( implicit_τy_coefficient(i, j, 1, grid, v_bottom_stress, clock, fields) / mᵢ * ℵᵢ 
        + implicit_τy_coefficient(i, j, 1, grid, v_top_stress, clock, fields) / mᵢ * ℵᵢ )

   Gⱽ = ifelse(mᵢ ≤ 0, zero(grid), Gⱽ)
   τᵢ = ifelse(mᵢ ≤ 0, zero(grid), τᵢ)

   return τᵢ, Gⱽ
end

@inline sum_of_forcing_x(i, j, 1, grid, rheology, u_forcing, fields, Δt) = u_forcing(i, j, 1, grid, fields)
@inline sum_of_forcing_y(i, j, 1, grid, rheology, u_forcing, fields, Δt) = v_forcing(i, j, 1, grid, fields)

# Fallback
@inline implicit_τx_coefficient(i, j, k, grid, stress, clock, fields) = zero(grid)
@inline implicit_τy_coefficient(i, j, k, grid, stress, clock, fields) = zero(grid)

# Fallback
@inline explicit_τx(i, j, k, grid, stress, clock, fields) = zero(grid)
@inline explicit_τy(i, j, k, grid, stress, clock, fields) = zero(grid)

@inline explicit_τx(i, j, k, grid, stress::Number, clock, fields) = stress
@inline explicit_τy(i, j, k, grid, stress::Number, clock, fields) = stress

@inline explicit_τx(i, j, k, grid, stress::AbstractArray, clock, fields) =  @inbounds stress[i, j, k] 
@inline explicit_τy(i, j, k, grid, stress::AbstractArray, clock, fields) =  @inbounds stress[i, j, k] 