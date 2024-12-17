using Oceananigans.Grids: AbstractGrid, architecture
using Oceananigans.Architectures: convert_args
using Oceananigans.Utils: configure_kernel
using Oceananigans.TimeSteppers: store_field_tendencies!
using Oceananigans.ImmersedBoundaries: retrieve_surface_active_cells_map

"""
    step_momentum!(model, rheology::AbstractExplicitRheology, Δt, χ)

function for stepping u and v in the case of _explicit_ solvers.
The sea-ice momentum equations are characterized by smaller time-scale than 
sea-ice thermodynamics and sea-ice tracer advection, therefore explicit rheologies require 
substepping over a set number of substeps.
"""
function step_momentum!(model, solver::ExplicitMomentumSolver, Δt, args...)

    grid = model.grid
    arch = architecture(grid)
    rheology = solver.rheology
    initialize_substepping!(model, solver)

    # The atmospheric stress component is fixed during time-stepping
    τua = model.external_momentum_stresses.u
    τva = model.external_momentum_stresses.v

    u, v = model.velocities

    immersed_bc = (; u = u.boundary_conditions.immersed, 
                     v = v.boundary_conditions.immersed)

    active_cells_map = retrieve_surface_active_cells_map(grid)

    # We step the momentum equation using a leap-frog scheme
    # where we alternate the order of solving u and v 
    args = (model.velocities, grid, 
            active_cells_map, Δt, 
            immersed_bc,
            model.clock,
            model.ocean_velocities,
            model.ocean_free_surface,
            model.coriolis,
            rheology,
            solver.auxiliary_fields,
            solver.substeps,
            model.ice_thickness,
            model.ice_concentration,
            model.ice_density,
            solver.ocean_ice_drag_coefficient,
            model.gravitational_acceleration)

    u_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _u_velocity_step!; active_cells_map)
    v_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _v_velocity_step!; active_cells_map)


    for substep in 1:solver.substeps
        Base.@_inline_meta
        # Compute stresses! depending on the particular rheology implementation
        compute_stresses!(model, solver, rheology, Δt)

        # TODO: This needs to be removed in some way!
        fill_rheology_halo_regions!(solver, rheology)

        # The momentum equations are solved using an alternating leap-frog algorithm
        # for u and v (used for the ocean - ice stresses and the coriolis term)
        # In even substeps we calculate uⁿ⁺¹ = f(vⁿ) and vⁿ⁺¹ = f(uⁿ⁺¹).
        # In odd substeps we switch and calculate vⁿ⁺¹ = f(uⁿ) and uⁿ⁺¹ = f(vⁿ⁺¹).
        if iseven(substep) 
            u_velocity_kernel!(args..., τua, nothing)
            v_velocity_kernel!(args..., τva, nothing)
        else
            v_velocity_kernel!(args..., τva, nothing)
            u_velocity_kernel!(args..., τua, nothing)
        end

        # TODO: This needs to be removed in some way!
        fill_halo_regions!(model.velocities)
    end

    return nothing
end
