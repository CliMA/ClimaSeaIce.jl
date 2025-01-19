using Oceananigans.Grids: AbstractGrid, architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Architectures: convert_args
using Oceananigans.Utils: configure_kernel
using Oceananigans.TimeSteppers: store_field_tendencies!
using Oceananigans.ImmersedBoundaries: retrieve_surface_active_cells_map, mask_immersed_field_xy!

struct SplitExplicitSolver 
    substeps :: Int
end

SplitExplicitSolver(; substeps=120) = SplitExplicitSolver(substeps)

const SplitExplicitMomentumEquation = SeaIceMomentumEquation{<:SplitExplicitSolver}

"""
    step_momentum!(model, rheology::AbstractExplicitRheology, Δt, χ)

function for stepping u and v in the case of _explicit_ solvers.
The sea-ice momentum equations are characterized by smaller time-scale than 
sea-ice thermodynamics and sea-ice tracer advection, therefore explicit rheologies require 
substepping over a set number of substeps.
"""
function step_momentum!(model, ice_dynamics::SplitExplicitMomentumEquation, Δt, args...)

    grid = model.grid
    arch = architecture(grid)
    rheology = ice_dynamics.rheology

    u, v = model.velocities
  
    args = (model.clock, 
            model.velocities, 
            ice_dynamics.coriolis, 
            model.ice_thickness, 
            model.ice_concentration, 
            model.ice_density)

    u_top_stress = model.external_momentum_stresses.top.u
    v_top_stress = model.external_momentum_stresses.top.v

    u_bottom_stress = model.external_momentum_stresses.bottom.u
    v_bottom_stress = model.external_momentum_stresses.bottom.v

    u_forcing = model.forcing.u
    v_forcing = model.forcing.v

    u_args = (args..., u_top_stress, u_bottom_stress, u_forcing)
    v_args = (args..., v_top_stress, v_bottom_stress, v_forcing)

    auxiliary_fields = ice_dynamics.auxiliary_fields

    u_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _u_velocity_step!)
    v_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _v_velocity_step!)

    substeps = ice_dynamics.solver.substeps
    
    fill_halo_regions!(model.velocities)
    initialize_rheology!(model, ice_dynamics.rheology)

    for substep in 1 : substeps
        # Compute stresses! depending on the particular rheology implementation
        compute_stresses!(model, ice_dynamics, rheology, Δt)

        # The momentum equations are solved using an alternating leap-frog algorithm
        # for u and v (used for the ocean - ice stresses and the coriolis term)
        # In even substeps we calculate uⁿ⁺¹ = f(vⁿ) and vⁿ⁺¹ = f(uⁿ⁺¹).
        # In odd substeps we switch and calculate vⁿ⁺¹ = f(uⁿ) and uⁿ⁺¹ = f(vⁿ⁺¹).
        if iseven(substep) 
            u_velocity_kernel!(u, grid, Δt, substeps, rheology, auxiliary_fields, u_args)
            v_velocity_kernel!(v, grid, Δt, substeps, rheology, auxiliary_fields, v_args)
        else
            v_velocity_kernel!(v, grid, Δt, substeps, rheology, auxiliary_fields, v_args)
            u_velocity_kernel!(u, grid, Δt, substeps, rheology, auxiliary_fields, u_args)
        end

        # TODO: This needs to be removed in some way!
        fill_halo_regions!(model.velocities)

        mask_immersed_field_xy!(model.velocities.u, k=1)
        mask_immersed_field_xy!(model.velocities.v, k=1)
    end

    return nothing
end

@kernel function _u_velocity_step!(u, grid, Δt, substeps, rheology, auxiliary_fields, args)
    i, j = @index(Global, NTuple)

    Δτ      = compute_time_stepᶠᶜᶜ(i, j, grid, Δt, rheology, substeps, auxiliary_fields) 
    τuᵢ, Gu = u_velocity_tendency(i, j, grid, Δτ, rheology, auxiliary_fields, args...)

    @inbounds u[i, j, 1] = (u[i, j, 1] + Δτ * Gu) / (1 + Δτ * τuᵢ)
end

@kernel function _v_velocity_step!(v, grid, Δt, substeps, rheology, auxiliary_fields, args)
    i, j = @index(Global, NTuple)
    
    Δτ      = compute_time_stepᶜᶠᶜ(i, j, grid, Δt, rheology, substeps, auxiliary_fields) 
    τvᵢ, Gv = v_velocity_tendency(i, j, grid, Δτ, rheology, auxiliary_fields, args...)
    
    @inbounds v[i, j, 1] = (v[i, j, 1] + Δτ * Gv) / (1 + Δτ * τvᵢ)
end