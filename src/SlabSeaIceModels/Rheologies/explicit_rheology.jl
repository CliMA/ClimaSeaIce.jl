using Oceananigans.Grids: AbstractGrid
using Oceananigans.TimeSteppers: store_field_tendencies!
using Printf

abstract type AbstractExplicitRheology <: AbstractRheology end

"""
    step_momentum!(model, rheology::AbstractExplicitRheology, Δt, χ)

function for stepping u and v in the case of _explicit_ rheologies.
The sea-ice momentum equations are characterized by smaller time-scale than 
sea-ice thermodynamics, therefore explicit rheologies require substepping
over a set number of substeps.
"""
function step_momentum!(model, rheology::AbstractExplicitRheology, Δt, χ)

    grid = model.grid
    arch = architecture(grid)

    Δτ = Δt / rheology.substeps

    τua = model.external_momentum_stress.u
    τva = model.external_momentum_stress.v

    # We step the momentum equation using a leap-frog scheme
    # where we alternate the order of solving u and v 
    for substep in 1:rheology.substeps
        
        compute_stresses!(model, model.rheology, Δτ)

        args = (model.velocities, grid, Δτ, 
                model.clock,
                model.ocean_velocities,
                model.coriolis,
                model.rheology,
                model.thickness)
                
        if iseven(substep)
            launch!(arch, grid, :xyz, _u_velocity_step!, args..., τua, nothing, fields(model))
            launch!(arch, grid, :xyz, _v_velocity_step!, args..., τva, nothing, fields(model))
        else
            launch!(arch, grid, :xyz, _v_velocity_step!, args..., τva, nothing, fields(model))
            launch!(arch, grid, :xyz, _u_velocity_step!, args..., τua, nothing, fields(model))
        end

        @info @sprintf
    end

    return nothing
end
