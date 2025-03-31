using Oceananigans.Grids: AbstractGrid, architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Utils: configure_kernel
using Oceananigans.ImmersedBoundaries: mask_immersed_field_xy!

struct SplitExplicitSolver 
    substeps :: Int
end

"""
    SplitExplicitSolver(; substeps=120)

Creates a `SplitExplicitSolver` that controls the dynamical evolution of sea-ice momentum
by subcycling `substeps` times in between each ice_thermodynamics / tracer advection time step.

The default number of substeps is 120.
"""
SplitExplicitSolver(; substeps=120) = SplitExplicitSolver(substeps)

const SplitExplicitMomentumEquation = SeaIceMomentumEquation{<:SplitExplicitSolver}

"""
    time_step_momentum!(model, rheology::AbstractExplicitRheology, Δt)

function for stepping u and v in the case of _explicit_ solvers.
The sea-ice momentum equations are characterized by smaller time-scale than 
sea-ice ice_thermodynamics and sea-ice tracer advection, therefore explicit rheologies require 
substepping over a set number of substeps.
"""
function time_step_momentum!(model, dynamics::SplitExplicitMomentumEquation, Δt)

    grid = model.grid
    arch = architecture(grid)
    rheology = dynamics.rheology

    u, v = model.velocities
  
    ocean_velocities = dynamics.ocean_velocities
    clock = model.clock
    coriolis = dynamics.coriolis

    minimum_mass = dynamics.minimum_mass
    minimum_concentration = dynamics.minimum_concentration

    top_stress = dynamics.external_momentum_stresses.top
    bottom_stress = dynamics.external_momentum_stresses.bottom

    u_forcing = model.forcing.u
    v_forcing = model.forcing.v
    u_immersed_bc = u.boundary_conditions.immersed
    v_immersed_bc = v.boundary_conditions.immersed

    model_fields = merge(dynamics.auxiliary_fields, model.velocities, 
                      (; h = model.ice_thickness, 
                         ℵ = model.ice_concentration, 
                         ρ = model.ice_density))

    u_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _u_velocity_step!)
    v_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _v_velocity_step!)

    substeps = dynamics.solver.substeps
    
    fill_halo_regions!(model.velocities)
    initialize_rheology!(model, dynamics.rheology)

    for substep in 1 : substeps
        # Compute stresses! depending on the particular rheology implementation
        compute_stresses!(model, dynamics, rheology, Δt)

        # The momentum equations are solved using an alternating leap-frog algorithm
        # for u and v (used for the ocean - ice stresses and the coriolis term)
        # In even substeps we calculate uⁿ⁺¹ = f(vⁿ) and vⁿ⁺¹ = f(uⁿ⁺¹).
        # In odd substeps we switch and calculate vⁿ⁺¹ = f(uⁿ) and uⁿ⁺¹ = f(vⁿ⁺¹).
        if iseven(substep) 
            u_velocity_kernel!(u, grid, Δt, substeps, rheology, model_fields, 
                               ocean_velocities, clock, coriolis,
                               minimum_mass, minimum_concentration, 
                               u_immersed_bc, top_stress, bottom_stress, u_forcing)

            v_velocity_kernel!(v, grid, Δt, substeps, rheology, model_fields, 
                               ocean_velocities, clock, coriolis, 
                               minimum_mass, minimum_concentration,
                               v_immersed_bc, top_stress, bottom_stress, v_forcing)

        else
            v_velocity_kernel!(v, grid, Δt, substeps, rheology, model_fields, 
                               ocean_velocities, clock, coriolis, 
                               minimum_mass, minimum_concentration,
                               v_immersed_bc, top_stress, bottom_stress, v_forcing)
            

            u_velocity_kernel!(u, grid, Δt, substeps, rheology, model_fields, 
                               ocean_velocities, clock, coriolis,
                               minimum_mass, minimum_concentration, 
                               u_immersed_bc, top_stress, bottom_stress, u_forcing)
        end

        # TODO: This needs to be removed in some way!
        fill_halo_regions!(model.velocities)

        mask_immersed_field_xy!(model.velocities.u, k=size(grid, 3))
        mask_immersed_field_xy!(model.velocities.v, k=size(grid, 3))
    end

    return nothing
end

@kernel function _u_velocity_step!(u, grid, Δt, 
                                   substeps, rheology, 
                                   model_fields, ocean_velocities, 
                                   clock, coriolis, 
                                   minimum_mass, minimum_concentration,
                                   u_immersed_bc, u_top_stress, u_bottom_stress, u_forcing)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) 

    mᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, ice_mass, model_fields.h, model_fields.ℵ, model_fields.ρ)
    ℵᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, model_fields.ℵ)

    Δτ = compute_substep_Δtᶠᶜᶜ(i, j, grid, Δt, rheology, substeps, model_fields) 
    Gu = u_velocity_tendency(i, j, grid, Δτ, rheology, model_fields, clock, coriolis, u_immersed_bc, u_top_stress, u_bottom_stress, u_forcing)
   
    # Implicit part of the stress that depends linearly on the velocity
    τuᵢ = ( implicit_τx_coefficient(i, j, kᴺ, grid, u_bottom_stress, clock, model_fields) 
          - implicit_τx_coefficient(i, j, kᴺ, grid, u_top_stress, clock, model_fields)) / mᵢ * ℵᵢ 

    τuᵢ = ifelse(mᵢ ≤ 0, zero(grid), τuᵢ)
    uᴰ  = @inbounds (u[i, j, 1] + Δτ * Gu) / (1 + Δτ * τuᵢ) # dynamical velocity 
    uᶠ  = free_drift_u(i, j, kᴺ, grid, ocean_velocities) # free drift velocity

    # If the ice mass or the ice concentration are below a certain threshold, 
    # the sea ice velocity is set to the free drift velocity
    sea_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)

    @inbounds u[i, j, 1] = ifelse(sea_ice, uᴰ, uᶠ)
end

@kernel function _v_velocity_step!(v, grid, Δt, 
                                   substeps, rheology, 
                                   model_fields, ocean_velocities, 
                                   clock, coriolis, 
                                   minimum_mass, minimum_concentration,
                                   v_immersed_bc, v_top_stress, v_bottom_stress, v_forcing)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) 

    mᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, ice_mass, model_fields.h, model_fields.ℵ, model_fields.ρ)
    ℵᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, model_fields.ℵ)
    
    Δτ = compute_substep_Δtᶜᶠᶜ(i, j, grid, Δt, rheology, substeps, model_fields) 
    Gv = v_velocity_tendency(i, j, grid, Δτ, rheology, model_fields, clock, coriolis, v_immersed_bc, v_top_stress, v_bottom_stress, v_forcing)

    # Implicit part of the stress that depends linearly on the velocity
    τvᵢ = ( implicit_τy_coefficient(i, j, kᴺ, grid, v_bottom_stress, clock, model_fields)
          - implicit_τy_coefficient(i, j, kᴺ, grid, v_top_stress, clock, model_fields)) / mᵢ * ℵᵢ 

    τvᵢ = ifelse(mᵢ ≤ 0, zero(grid), τvᵢ)

    vᴰ = @inbounds (v[i, j, 1] + Δτ * Gv) / (1 + Δτ * τvᵢ)# dynamical velocity 
    vᶠ = free_drift_v(i, j, kᴺ, grid, ocean_velocities) # free drift velocity

    # If the ice mass or the ice concentration are below a certain threshold, 
    # the sea ice velocity is set to the free drift velocity
    sea_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)

    @inbounds v[i, j, 1] = ifelse(sea_ice, vᴰ, vᶠ)
end
