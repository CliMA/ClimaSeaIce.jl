using Oceananigans.Grids: AbstractGrid, architecture
using Oceananigans.Architectures: convert_args
using Oceananigans.Utils: configure_kernel
using Oceananigans.TimeSteppers: store_field_tendencies!
using Oceananigans.ImmersedBoundaries: retrieve_surface_active_cells_map

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
            ice_dynamics.rheology, 
            ice_dynamics.auxiliary_fields, 
            model.ice_thickness, 
            model.ice_concentration, 
            model.ice_density)

    u_top_stress = model.external_momentum_stresses.top.u
    v_top_stress = model.external_momentum_stresses.top.v
    u_bottom_stress = model.external_momentum_stresses.bottom.u
    v_bottom_stress = model.external_momentum_stresses.bottom.v

    u_args = (args..., u_top_stress, u_bottom_stress)
    v_args = (args..., v_top_stress, v_bottom_stress)

    u_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _u_velocity_step!)
    v_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _v_velocity_step!)

    for substep in 1:ice_dynamics.solver.substeps
        # Compute stresses! depending on the particular rheology implementation
        compute_stresses!(model, ice_dynamics, rheology, Δt)

        # The momentum equations are solved using an alternating leap-frog algorithm
        # for u and v (used for the ocean - ice stresses and the coriolis term)
        # In even substeps we calculate uⁿ⁺¹ = f(vⁿ) and vⁿ⁺¹ = f(uⁿ⁺¹).
        # In odd substeps we switch and calculate vⁿ⁺¹ = f(uⁿ) and uⁿ⁺¹ = f(vⁿ⁺¹).
        if iseven(substep) 
            u_velocity_kernel!(u, grid, u_args)
            v_velocity_kernel!(v, grid, v_args)
        else
            v_velocity_kernel!(v, grid, v_args)
            u_velocity_kernel!(u, grid, u_args)
        end

        # TODO: This needs to be removed in some way!
        fill_halo_regions!(model.velocities)
    end

    return nothing
end

@kernel function _u_velocity_step!(u, grid, args)
    i, j = @index(Global, NTuple)
    Gⁿu = u_velocity_tendency(i, j, grid, args...)
    @inbounds u[i, j, 1] += Δt * Gⁿu
end

@kernel function _v_velocity_step!(v, grid, args)
    i, j = @index(Global, NTuple)
    Gⁿv = v_velocity_tendency(i, j, grid, args...)
    @inbounds v[i, j, 1] += Δt * Gⁿv
end