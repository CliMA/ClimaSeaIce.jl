using Oceananigans.Utils
using ClimaSeaIce: timestepping_coefficients

const ExplicitMomentumEquation = SeaIceMomentumEquation{<:ExplicitSolver}

# Simple explicit stepping of the momentum equations
function step_momentum!(model, ::ExplicitMomentumEquation, Δt, stage)
    grid = model.grid
    arch = architecture(grid)

    u, v = model.velocities
    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻

    dynamics = model.dynamics

    α, β = timestepping_coefficients(model.timestepper, stage)

    model_fields = merge(dynamics.auxiliary_fields, model.velocities, 
                      (; h = model.ice_thickness, 
                         ℵ = model.ice_concentration, 
                         ρ = model.ice_density))


    ocean_velocities = dynamics.ocean_velocities
    clock = model.clock
    minimum_mass = dynamics.minimum_mass
    minimum_concentration = dynamics.minimum_concentration

    top_stress = dynamics.external_momentum_stresses.top
    bottom_stress = dynamics.external_momentum_stresses.bottom

    launch!(arch, grid, :xyz, _step_velocities!, u, v, grid, Gⁿ, G⁻, Δt, α, β,
            τ_top, τ_bottom, ocean_velocities, 
            minimum_mass, minimum_concentration, clock, model_fields)

    return nothing
end

@kernel function _step_velocities!(u, v, grid, Gⁿ, G⁻, Δt, α, β, 
                                   top_stress, bottom_stress, 
                                   ocean_velocities, minimum_mass, minimum_concentration, clock, fields)

    i, j, k = @index(Global, NTuple)

    ℵᶠᶜ = ℑxᶠᵃᵃ(i, j, k, grid, fields.ℵ)
    ℵᶜᶠ = ℑyᵃᶠᵃ(i, j, k, grid, fields.ℵ)
    mᶠᶜ = ℑxᶠᵃᵃ(i, j, k, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)
    mᶜᶠ = ℑyᵃᶠᵃ(i, j, k, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)

   # Implicit part of the stress that depends linearly on the velocity
   τuᵢ = ( implicit_τx_coefficient(i, j, 1, grid, bottom_stress, clock, fields) 
         + implicit_τx_coefficient(i, j, 1, grid, top_stress,    clock, fields)) / mᶠᶜ * ℵᶠᶜ 
   
   τvᵢ = ( implicit_τy_coefficient(i, j, 1, grid, bottom_stress, clock, fields) 
         + implicit_τy_coefficient(i, j, 1, grid, top_stress,    clock, fields)) / mᶜᶠ * ℵᶜᶠ 

    @inbounds begin
        uᴰ = (Δt * (α * Gⁿ.u[i, j, k] + β * G⁻.u[i, j, k])) / (1 + Δt * τuᵢ)
        vᴰ = (Δt * (α * Gⁿ.v[i, j, k] + β * G⁻.v[i, j, k])) / (1 + Δt * τvᵢ)

        uᶠ = free_drift_u(i, j, k, grid, ocean_velocities)
        vᶠ = free_drift_v(i, j, k, grid, ocean_velocities)

        sea_ice = (mᶠᶜ ≥ minimum_mass) & (ℵᶠᶜ ≥ minimum_concentration)
        u[i, j, k] = ifelse(sea_ice, uᴰ, uᶠ)

        sea_ice = (mᶜᶠ ≥ minimum_mass) & (ℵᶜᶠ ≥ minimum_concentration)
        v[i, j, k] = ifelse(sea_ice, vᴰ, vᶠ)
    end 
end

# Compute the tendencies for the explicit momentum equations
function compute_momentum_tendencies!(model, ::ExplicitMomentumEquation, Δt)
    
    dynamics = model.dynamics
    grid = model.grid

    clock    = model.clock
    coriolis = dynamics.coriolis
    rheology = dynamics.rheology

    model_fields = merge(dynamics.auxiliary_fields, model.velocities, 
            (; h = model.ice_thickness, 
               ℵ = model.ice_concentration, 
               ρ = model.ice_density))

    u_top_stress = dynamics.external_momentum_stresses.top.u
    v_top_stress = dynamics.external_momentum_stresses.top.v
    u_bottom_stress = dynamics.external_momentum_stresses.bottom.u
    v_bottom_stress = dynamics.external_momentum_stresses.bottom.v

    Gu = model.timestepper.Gⁿ.u
    Gv = model.timestepper.Gⁿ.v

    launch!(architecture(grid), grid, :xy, _compute_velocity_tendencies!, Gu, Gv, grid, Δt,
            rheology, model_fields, clock, coriolis,
            u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress, model.forcing)

    return nothing
end

@kernel function _compute_velocity_tendencies!(Gu, Gv, grid, Δt,
                                               rheology, model_fields, clock, coriolis,
                                               u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress, forcing)
    i, j = @index(Global, NTuple)
    @inbounds Gu[i, j, 1] = u_velocity_tendency(i, j, grid, Δt,
                                                rheology, model_fields, clock, coriolis,
                                                u_top_stress, u_bottom_stress, forcing.u)

    @inbounds Gv[i, j, 1] = v_velocity_tendency(i, j, grid, Δt,
                                                rheology, model_fields, clock, coriolis,
                                                v_top_stress, v_bottom_stress, forcing.v)
end