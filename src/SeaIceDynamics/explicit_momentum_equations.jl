using Oceananigans.Utils

const ExplicitMomentumEquation = SeaIceMomentumEquation{<:ExplicitSolver}

# Simple explicit stepping of the momentum equations
function time_step_momentum!(model, ::ExplicitMomentumEquation, Δt)
    grid = model.grid
    arch = architecture(grid)

    u, v = model.velocities
    Gⁿ = model.timestepper.Gⁿ

    dynamics = model.dynamics

    model_fields = merge(dynamics.auxiliary_fields, model.velocities, 
                      (; h = model.ice_thickness, 
                         ℵ = model.ice_concentration, 
                         ρ = model.ice_density))


    free_drift = dynamics.free_drift
    clock = model.clock
    minimum_mass = dynamics.minimum_mass
    minimum_concentration = dynamics.minimum_concentration

    top_stress = dynamics.external_momentum_stresses.top
    bottom_stress = dynamics.external_momentum_stresses.bottom

    launch!(arch, grid, :xy, _step_velocities!, u, v, grid, Gⁿ, Δt, 
            top_stress, bottom_stress, free_drift, 
            minimum_mass, minimum_concentration, clock, model_fields)

    return nothing
end

@kernel function _step_velocities!(u, v, grid, Gⁿ, Δt, 
                                   top_stress, bottom_stress, 
                                   free_drift, minimum_mass, minimum_concentration, clock, fields)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    ℵᶠᶠ = ℑxyᶠᶠᵃ(i, j, 1, grid, fields.ℵ)
    mᶠᶠ = ℑxyᶠᶠᵃ(i, j, 1, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)
    
   # Implicit part of the stress that depends linearly on the velocity
   τuᵢ = ( implicit_τx_coefficient(i, j, kᴺ, grid, bottom_stress, clock, fields) 
         - implicit_τx_coefficient(i, j, kᴺ, grid, top_stress,    clock, fields)) / mᶠᶠ * ℵᶠᶠ 
   
   τvᵢ = ( implicit_τy_coefficient(i, j, kᴺ, grid, bottom_stress, clock, fields) 
         - implicit_τy_coefficient(i, j, kᴺ, grid, top_stress,    clock, fields)) / mᶠᶠ * ℵᶠᶠ 

    @inbounds begin
        uᴰ = (u[i, j, 1] + Δt * Gⁿ.u[i, j, 1]) / (1 + Δt * τuᵢ)
        vᴰ = (v[i, j, 1] + Δt * Gⁿ.v[i, j, 1]) / (1 + Δt * τvᵢ)

        uᶠ = free_drift_u(i, j, kᴺ, grid, free_drift, clock, fields)
        vᶠ = free_drift_v(i, j, kᴺ, grid, free_drift, clock, fields)

        sea_ice = (mᶠᶠ ≥ minimum_mass) & (ℵᶠᶠ ≥ minimum_concentration)
        
        u[i, j, 1] = ifelse(sea_ice, uᴰ, uᶠ)
        v[i, j, 1] = ifelse(sea_ice, vᴰ, vᶠ)
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

    top_stress = dynamics.external_momentum_stresses.top
    bottom_stress = dynamics.external_momentum_stresses.bottom

    Gu = model.timestepper.Gⁿ.u
    Gv = model.timestepper.Gⁿ.v

    launch!(architecture(grid), grid, :xy, _compute_velocity_tendencies!, Gu, Gv, grid, Δt,
            rheology, model_fields, clock, coriolis,
            top_stress, bottom_stress, model.forcing)

    return nothing
end

@kernel function _compute_velocity_tendencies!(Gu, Gv, grid, Δt,
                                               rheology, model_fields, clock, coriolis,
                                               top_stress, bottom_stress, forcing)
    i, j = @index(Global, NTuple)
    @inbounds Gu[i, j, 1] = u_velocity_tendency(i, j, grid, Δt,
                                                rheology, model_fields, clock, coriolis,
                                                top_stress, bottom_stress, forcing.u)

    @inbounds Gv[i, j, 1] = v_velocity_tendency(i, j, grid, Δt,
                                                rheology, model_fields, clock, coriolis,
                                                top_stress, bottom_stress, forcing.v)
end