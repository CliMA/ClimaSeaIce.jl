using Oceananigans.Grids: AbstractGrid, architecture
using Oceananigans.Architectures: convert_args
using Oceananigans.TimeSteppers: store_field_tendencies!
using ClimaSeaIce.SeaIceDynamics: AbstractRheology
using Printf

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

    velocity_boundary_conditions = (; u = u.boundary_conditions, 
                                      v = v.boundary_conditions)

    # We step the momentum equation using a leap-frog scheme
    # where we alternate the order of solving u and v 
    args = (model.velocities, grid, Δt, 
            velocity_boundary_conditions,
            model.clock,
            model.ocean_velocities,
            model.coriolis,
            rheology,
            solver.auxiliary_fields,
            solver.substeps,
            solver.substepping_coefficient,
            model.ice_thickness,
            model.ice_concentration,
            model.ice_density,
            model.ocean_density,
            solver.ocean_ice_drag_coefficient)

    GC.@preserve args begin
        # We need to perform ~100 time-steps which means
        # launching ~200 very small kernels: we are limited by
        # latency of argument conversion to GPU-compatible values.
        # To alleviate this penalty we convert first and then we substep!
        converted_args = convert_args(arch, args)

        for substep in 1:solver.substeps
            Base.@_inline_meta
            # Compute stresses! depending on the particular rheology implementation
            compute_stresses!(model, solver, rheology, Δt)

            # The momentum equations are solved using an alternating leap-frog algorithm
            # for u and v (used for the ocean - ice stresses and the coriolis term)
            # In even substeps we calculate uⁿ⁺¹ = f(vⁿ) and vⁿ⁺¹ = f(uⁿ⁺¹).
            # In odd substeps we switch and calculate vⁿ⁺¹ = f(uⁿ) and uⁿ⁺¹ = f(vⁿ⁺¹).
            if iseven(substep) 
                launch!(arch, grid, :xy, _u_velocity_step!, converted_args..., τua, nothing, fields(model))
                launch!(arch, grid, :xy, _v_velocity_step!, converted_args..., τva, nothing, fields(model))
            else
                launch!(arch, grid, :xy, _v_velocity_step!, converted_args..., τva, nothing, fields(model))
                launch!(arch, grid, :xy, _u_velocity_step!, converted_args..., τua, nothing, fields(model))
            end
        end
    end

    return nothing
end
