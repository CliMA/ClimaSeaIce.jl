using Oceananigans.Grids: AbstractGrid, architecture
using Oceananigans.TimeSteppers: store_field_tendencies!
using ClimaSeaIce.SlabSeaIceModels.SlabSeaIceDynamics: AbstractRheology
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

    initialize_substepping!(model, rheology)

    # The atmospheric stress component is fixed during time-stepping
    τua = model.external_momentum_stress.u
    τva = model.external_momentum_stress.v

    # We step the momentum equation using a leap-frog scheme
    # where we alternate the order of solving u and v 
    for substep in 1:rheology.substeps
        
        # Compute stresses! depends on the particular 
        # rheology implemented
        compute_stresses!(model, model.rheology, Δt)

        args = (model.velocities, grid, Δt, 
                model.clock,
                model.ocean_velocities,
                model.coriolis,
                model.rheology,
                model.ice_thickness,
                model.concentration)
                
        # The momentum equations are solved using an alternating leap-frog algorithm
        # for u and v (used for the ocean - ice stresses and the coriolis term)
        # In even substeps we calculate uⁿ⁺¹ = f(vⁿ) and vⁿ⁺¹ = f(uⁿ⁺¹).
        # In odd substeps we swith and calculate vⁿ⁺¹ = f(uⁿ) and uⁿ⁺¹ = f(vⁿ⁺¹).
        if iseven(substep)
            launch!(arch, grid, :xyz, _u_velocity_step!, args..., τua, nothing, fields(model))
            launch!(arch, grid, :xyz, _v_velocity_step!, args..., τva, nothing, fields(model))
        else
            launch!(arch, grid, :xyz, _v_velocity_step!, args..., τva, nothing, fields(model))
            launch!(arch, grid, :xyz, _u_velocity_step!, args..., τua, nothing, fields(model))
        end
    end

    return nothing
end
