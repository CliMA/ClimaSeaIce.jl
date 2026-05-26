using Oceananigans.Utils

const ExplicitMomentumEquation = SeaIceMomentumEquation{<:ExplicitSolver}

previous_velocities(model, timestepper) = model.velocities
previous_velocities(model, timestepper::SplitRungeKuttaTimeStepper) = (u = model.timestepper.Ψ⁻.u, v = model.timestepper.Ψ⁻.v)

# Simple explicit stepping of the momentum equations
function time_step_momentum!(model, ::ExplicitMomentumEquation, Δt)
    grid = model.grid
    arch = architecture(grid)

    u,  v  = model.velocities
    u⁻, v⁻ = previous_velocities(model, model.timestepper)
    Gⁿ = model.timestepper.Gⁿ

    dynamics = model.dynamics

    model_fields = merge(dynamics.auxiliaries.fields, model.velocities, 
                      (; h = model.ice_thickness, 
                         ℵ = model.ice_concentration, 
                         ρ = model.sea_ice_density))
                
    free_drift = dynamics.free_drift
    clock = model.clock
    minimum_mass = dynamics.minimum_mass
    minimum_concentration = dynamics.minimum_concentration

    top_stress = dynamics.external_momentum_stresses.top
    bottom_stress = dynamics.external_momentum_stresses.bottom

    launch!(arch, grid, :xy, _step_velocities!, u, v, u⁻, v⁻, grid, Gⁿ, Δt, 
            top_stress, bottom_stress, free_drift, 
            minimum_mass, minimum_concentration, clock, model_fields)

    return nothing
end

@kernel function _step_velocities!(u, v, u⁻, v⁻, grid, Gⁿ, Δt, 
                                   top_stress, bottom_stress, 
                                   free_drift, minimum_mass, minimum_concentration, clock, fields)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)
    ℵᶠᶜ  = ℑxᶠᵃᵃ(i, j, kᴺ, grid, fields.ℵ)
    mᶠᶜ  = ℑxᶠᵃᵃ(i, j, kᴺ, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)
    ℵᶜᶠ  = ℑyᵃᶠᵃ(i, j, kᴺ, grid, fields.ℵ)
    mᶜᶠ  = ℑyᵃᶠᵃ(i, j, kᴺ, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)

   # Implicit part of the stress that depends linearly on the velocity. Read the value
   # precomputed by `compute_implicit_stress_coefficients!`, since this kernel updates u and v
   # together and an on-the-fly evaluation would race on the neighbour velocities.
   τuᵢ = ( implicit_τx_coefficient_field(i, j, kᴺ, grid, bottom_stress, clock, fields)
         - implicit_τx_coefficient_field(i, j, kᴺ, grid, top_stress,    clock, fields)) / mᶠᶜ * ℵᶠᶜ

   τvᵢ = ( implicit_τy_coefficient_field(i, j, kᴺ, grid, bottom_stress, clock, fields)
         - implicit_τy_coefficient_field(i, j, kᴺ, grid, top_stress,    clock, fields)) / mᶜᶠ * ℵᶜᶠ

    @inbounds begin
        uᴰ = (u⁻[i, j, 1] + Δt * Gⁿ.u[i, j, 1]) / (1 + Δt * τuᵢ)
        uᶠ = free_drift_u(i, j, kᴺ, grid, free_drift, clock, fields)
        
        sea_ice = (mᶠᶜ ≥ minimum_mass) & (ℵᶠᶜ ≥ minimum_concentration)
        u[i, j, 1] = ifelse(sea_ice, uᴰ, uᶠ)

        vᴰ = (v⁻[i, j, 1] + Δt * Gⁿ.v[i, j, 1]) / (1 + Δt * τvᵢ)
        vᶠ = free_drift_v(i, j, kᴺ, grid, free_drift, clock, fields)
       
        sea_ice = (mᶜᶠ ≥ minimum_mass) & (ℵᶜᶠ ≥ minimum_concentration)
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

    model_fields = merge(dynamics.auxiliaries.fields, model.velocities, 
            (; h = model.ice_thickness, 
               ℵ = model.ice_concentration, 
               ρ = model.sea_ice_density))

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
