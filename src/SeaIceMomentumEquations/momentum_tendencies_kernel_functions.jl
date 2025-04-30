using Oceananigans.Coriolis: y_f_cross_U, x_f_cross_U

"""compute explicit ice u-velocity tendencies"""
@inline function u_velocity_tendency(i, j, grid, Δt,
                                     rheology,
                                     model_fields,
                                     clock,
                                     coriolis,
                                     u_immersed_bc,
                                     u_top_stress,
                                     u_bottom_stress,
                                     u_forcing)

     kᴺ = size(grid, 3)
     h  = model_fields.h
     ℵ  = model_fields.ℵ
     ρ  = model_fields.ρ
     U  = (u = model_fields.u, v = model_fields.v)

     # Ice mass (per unit area) interpolated on u points
     ℵᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, ℵ)
     mᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, ice_mass, h, ℵ, ρ) 

     Gᵁ = ( - x_f_cross_U(i, j, kᴺ, grid, coriolis, U) 
            - explicit_τx(i, j, kᴺ, grid, u_top_stress, clock, model_fields) / mᵢ * ℵᵢ
            + explicit_τx(i, j, kᴺ, grid, u_bottom_stress, clock, model_fields) / mᵢ * ℵᵢ
            + ∂ⱼ_σ₁ⱼ(i, j, kᴺ, grid, rheology, clock, model_fields) / mᵢ
            + immersed_∂ⱼ_σ₁ⱼ(i, j, kᴺ, grid, u_immersed_bc, rheology, clock, model_fields) / mᵢ
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
                                     v_immersed_bc,
                                     v_top_stress,
                                     v_bottom_stress,
                                     v_forcing)

     kᴺ = size(grid, 3)
     h = model_fields.h
     ℵ = model_fields.ℵ
     ρ = model_fields.ρ
     U = (u = model_fields.u, v = model_fields.v)

     # Ice mass (per unit area) interpolated on v points
     ℵᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, ℵ)
     mᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, ice_mass, h, ℵ, ρ) 

     Gⱽ = ( - y_f_cross_U(i, j, kᴺ, grid, coriolis, U)
            - explicit_τy(i, j, kᴺ, grid, v_top_stress, clock, model_fields) / mᵢ * ℵᵢ
            + explicit_τy(i, j, kᴺ, grid, v_bottom_stress, clock, model_fields) / mᵢ * ℵᵢ
            + ∂ⱼ_σ₂ⱼ(i, j, kᴺ, grid, rheology, clock, model_fields) / mᵢ 
            + immersed_∂ⱼ_σ₂ⱼ(i, j, kᴺ, grid, v_immersed_bc, rheology, clock, model_fields) / mᵢ
            + sum_of_forcing_v(i, j, kᴺ, grid, rheology, v_forcing, model_fields, Δt)) # sum of user defined forcing and possibly other forcing terms that are rheology-dependent 

     Gⱽ = ifelse(mᵢ ≤ 0, zero(grid), Gⱽ)

     return Gⱽ
end
