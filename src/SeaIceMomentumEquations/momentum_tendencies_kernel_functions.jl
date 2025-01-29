using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U
using Oceananigans.ImmersedBoundaries: active_linear_index_to_tuple

"""compute ice u-velocity tendencies"""
@inline function u_velocity_tendency(i, j, grid, Δt,
                                     rheology,
                                     model_fields,
                                     clock,
                                     coriolis,
                                     u_top_stress,
                                     u_bottom_stress,
                                     u_forcing)

     h = model_fields.h
     ℵ = model_fields.ℵ
     ρ = model_fields.ρ

     # Ice mass (per unit area) interpolated on u points
     ℵᵢ = ℑxᶠᵃᵃ(i, j, 1, grid, ℵ)
     mᵢ = ℑxᶠᵃᵃ(i, j, 1, grid, ice_mass, h, ρ) 

     Gᵁ = ( - x_f_cross_U(i, j, 1, grid, coriolis, model_fields) 
            + explicit_τx(i, j, 1, grid, u_top_stress, clock, model_fields) / mᵢ * ℵᵢ
            + explicit_τx(i, j, 1, grid, u_bottom_stress, clock, model_fields) / mᵢ * ℵᵢ
            + ∂ⱼ_σ₁ⱼ(i, j, 1, grid, rheology, clock, model_fields, Δt) / mᵢ
            # sum of user defined forcing and possibly other forcing terms that are rheology-dependent 
            + sum_of_forcing_x(i, j, 1, grid, rheology, u_forcing, model_fields, Δt)) 

     # Implicit part of the stress that depends linearly on the velocity
     τᵢ = ( implicit_τx_coefficient(i, j, 1, grid, u_bottom_stress, clock, model_fields) / mᵢ * ℵᵢ
          + implicit_τx_coefficient(i, j, 1, grid, u_top_stress, clock, model_fields) / mᵢ * ℵᵢ )

     Gᵁ = ifelse(mᵢ ≤ 0, zero(grid), Gᵁ)
     τᵢ = ifelse(mᵢ ≤ 0, zero(grid), τᵢ)

     return τᵢ, Gᵁ
end

"""compute ice v-velocity tendencies"""
@inline function v_velocity_tendency(i, j, grid, Δt,
                                     rheology,
                                     model_fields,
                                     clock,
                                     coriolis,
                                     v_top_stress,
                                     v_bottom_stress,
                                     v_forcing)

     h = model_fields.h
     ℵ = model_fields.ℵ
     ρ = model_fields.ρ

     # Ice mass (per unit area) interpolated on v points
     ℵᵢ = ℑyᵃᶠᵃ(i, j, 1, grid, ℵ)
     mᵢ = ℑyᵃᶠᵃ(i, j, 1, grid, ice_mass, h, ρ) 

     Gⱽ = ( - y_f_cross_U(i, j, 1, grid, coriolis, model_fields)
            + explicit_τy(i, j, 1, grid, v_top_stress, clock, model_fields) / mᵢ * ℵᵢ
            + explicit_τy(i, j, 1, grid, v_bottom_stress, clock, model_fields) / mᵢ * ℵᵢ
            + ∂ⱼ_σ₂ⱼ(i, j, 1, grid, rheology, clock, model_fields, Δt) / mᵢ 
            # sum of user defined forcing and possibly other forcing terms that are rheology-dependent 
            + sum_of_forcing_y(i, j, 1, grid, rheology, v_forcing, model_fields, Δt))

     # Implicit part of the stress that depends linearly on the velocity
     τᵢ = ( implicit_τy_coefficient(i, j, 1, grid, v_bottom_stress, clock, model_fields) / mᵢ * ℵᵢ 
          + implicit_τy_coefficient(i, j, 1, grid, v_top_stress, clock, model_fields) / mᵢ * ℵᵢ )

     Gⱽ = ifelse(mᵢ ≤ 0, zero(grid), Gⱽ)
     τᵢ = ifelse(mᵢ ≤ 0, zero(grid), τᵢ)

   return τᵢ, Gⱽ
end

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