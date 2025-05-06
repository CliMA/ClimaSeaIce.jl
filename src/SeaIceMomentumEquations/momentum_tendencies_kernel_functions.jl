using Oceananigans.Coriolis: fᶠᶠᵃ

"""compute explicit ice u-velocity tendencies"""
@inline function u_velocity_tendency(i, j, grid, Δt,
                                     rheology,
                                     model_fields,
                                     clock,
                                     coriolis,
                                     u_top_stress,
                                     u_bottom_stress,
                                     u_forcing)

     kᴺ = size(grid, 3)
     h = model_fields.h
     ℵ = model_fields.ℵ
     ρ = model_fields.ρ

     # Ice mass (per unit area) interpolated on u points
     ℵᵢ = ℑxyᶠᶠᵃ(i, j, 1, grid, ℵ)
     mᵢ = ℑxyᶠᶠᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρ) 

     @inbounds Gᵁ = ( + fᶠᶠᵃ(i, j, 1, grid, coriolis) * model_fields.v[i, j, 1]
                      - explicit_τx(i, j, kᴺ, grid, u_top_stress, clock, model_fields) / mᵢ * ℵᵢ
                      + explicit_τx(i, j, kᴺ, grid, u_bottom_stress, clock, model_fields) / mᵢ * ℵᵢ
                      + ∂ⱼ_σ₁ⱼ(i, j, 1, grid, rheology, clock, model_fields) / mᵢ
                      + sum_of_forcing_u(i, j, kᴺ, grid, rheology, u_forcing, model_fields, Δt))  # sum of user defined forcing and possibly other forcing terms that are rheology-dependent 

     Gᵁ = ifelse(mᵢ ≤ 0, zero(grid), Gᵁ)

     return Gᵁ
end

"""compute explicit ice v-velocity tendencies"""
@inline function v_velocity_tendency(i, j, grid, Δt,
                                     rheology,
                                     model_fields,
                                     clock,
                                     coriolis,
                                     v_top_stress,
                                     v_bottom_stress,
                                     v_forcing)

     kᴺ = size(grid, 3)
     h = model_fields.h
     ℵ = model_fields.ℵ
     ρ = model_fields.ρ

     # Ice mass (per unit area) interpolated on v points
     ℵᵢ = ℑxyᶠᶠᵃ(i, j, 1, grid, ℵ)
     mᵢ = ℑxyᶠᶠᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρ) 

     @inbounds Gⱽ = ( - fᶠᶠᵃ(i, j, 1, grid, coriolis) * model_fields.u[i, j, 1]
                      - explicit_τy(i, j, kᴺ, grid, v_top_stress, clock, model_fields) / mᵢ * ℵᵢ
                      + explicit_τy(i, j, kᴺ, grid, v_bottom_stress, clock, model_fields) / mᵢ * ℵᵢ
                      + ∂ⱼ_σ₂ⱼ(i, j, 1, grid, rheology, clock, model_fields) / mᵢ 
                      + sum_of_forcing_v(i, j, kᴺ, grid, rheology, v_forcing, model_fields, Δt)) # sum of user defined forcing and possibly other forcing terms that are rheology-dependent 

     Gⱽ = ifelse(mᵢ ≤ 0, zero(grid), Gⱽ)

     return Gⱽ
end
