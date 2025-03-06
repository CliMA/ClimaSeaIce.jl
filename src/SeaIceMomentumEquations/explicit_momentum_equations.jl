using Oceananigans.Utils

const ExplicitMomentumEquation = SeaIceMomentumEquation{<:ExplicitSolver}

# Simple explicit stepping of the momentum equations
function step_momentum!(model, ::ExplicitMomentumEquation, Δt, args...)
    grid = model.grid
    arch = architecture(grid)

    u, v = model.velocities
    Gⁿ = model.timestepper.Gⁿ

    dynamics = model.dynamics

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

    launch!(arch, grid, :xy, _step_velocities!, u, v, grid, Gⁿ, Δt, 
            top_stress, bottom_stress, ocean_velocities, 
            minimum_mass, minimum_concentration, clock, model_fields)

    return nothing
end

@kernel function _step_velocities!(u, v, grid, Gⁿ, Δt, 
                                   top_stress, bottom_stress, 
                                   ocean_velocities, minimum_mass, minimum_concentration, clock, fields)

    i, j = @index(Global, NTuple)
    
    ℵᶠᶠ = ℑxyᶠᶠᵃ(i, j, 1, grid, fields.ℵ)
    mᶠᶠ = ℑxyᶠᶠᵃ(i, j, 1, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)
    
   # Implicit part of the stress that depends linearly on the velocity
   τuᵢ = ( implicit_τx_coefficient(i, j, 1, grid, bottom_stress, clock, fields) 
         + implicit_τx_coefficient(i, j, 1, grid, top_stress,    clock, fields)) / mᶠᶜ * ℵᶠᶜ 
   
   τvᵢ = ( implicit_τy_coefficient(i, j, 1, grid, bottom_stress, clock, fields) 
         + implicit_τy_coefficient(i, j, 1, grid, top_stress,    clock, fields)) / mᶜᶠ * ℵᶜᶠ 

    @inbounds begin
        uᴰ = (u[i, j, 1] + Δt * Gⁿ.u[i, j, 1]) / (1 + Δt * τuᵢ)
        vᴰ = (v[i, j, 1] + Δt * Gⁿ.v[i, j, 1]) / (1 + Δt * τvᵢ)

        uᶠ = free_drift_u(i, j, 1, grid, ocean_velocities)
        vᶠ = free_drift_v(i, j, 1, grid, ocean_velocities)

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

    u_immersed_bc = model_fields.u.boundary_conditions.immersed
    v_immersed_bc = model_fields.v.boundary_conditions.immersed

    Gu = model.timestepper.Gⁿ.u
    Gv = model.timestepper.Gⁿ.v

    launch!(architecture(grid), grid, :xy, _compute_velocity_tendencies!, Gu, Gv, grid, Δt,
            rheology, model_fields, clock, coriolis,
            u_immersed_bc, v_immersed_bc,
            top_stress, bottom_stress, model.forcing)

    return nothing
end

@kernel function _compute_velocity_tendencies!(Gu, Gv, grid, Δt,
                                               rheology, model_fields, clock, coriolis,
                                               u_immersed_bc, v_immersed_bc,
                                               top_stress, bottom_stress, forcing)
    i, j = @index(Global, NTuple)
    @inbounds Gu[i, j, 1] = u_velocity_tendency(i, j, grid, Δt,
                                                rheology, model_fields, clock, coriolis,
                                                u_immersed_bc, top_stress, bottom_stress, forcing.u)

    @inbounds Gv[i, j, 1] = v_velocity_tendency(i, j, grid, Δt,
                                                rheology, model_fields, clock, coriolis,
                                                v_immersed_bc, top_stress, bottom_stress, forcing.v)
end