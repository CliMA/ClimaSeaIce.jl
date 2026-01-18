#TODO: Move to Oceananigans??
using Oceananigans.TimeSteppers: AbstractTimeStepper
using Oceananigans.Fields: FunctionField, location
using Oceananigans.Utils: @apply_regionally, apply_regionally!

import Oceananigans: prognostic_state, restore_prognostic_state!

mutable struct ForwardEulerTimeStepper{FT, GT, IT} <: AbstractTimeStepper
                 Gⁿ :: GT
    implicit_solver :: IT
end

"""
        ForwardEulerTimeStepper(grid, prognostic_fields;
                                implicit_solver = nothing,
                                Gⁿ = map(similar, prognostic_fields))

Return a 1st-order Forward-Euler (FE) time stepper (`ForwardEulerTimeStepper`)
on `grid`, with `tracers`. The tendency fields `Gⁿ`, usually equal to
the prognostic_fields passed as positional argument, can be specified via  optional `kwargs`.

The 1st-order Forward-Euler timestepper steps forward the state `Uⁿ` by `Δt` via

```julia
Uⁿ⁺¹ = Uⁿ + Δt * Gⁿ

where `Uⁿ` is the state at the ``n``-th timestep and `Gⁿ` is the tendency
at the ``n``-th timestep.
"""
function ForwardEulerTimeStepper(grid, prognostic_fields;
                                 implicit_solver::IT = nothing,
                                 Gⁿ = map(similar, prognostic_fields)) where IT

    FT = eltype(grid)
    GT = typeof(Gⁿ)

    return ForwardEulerTimeStepper{FT, GT, IT}(Gⁿ, implicit_solver)
end

reset!(timestepper::ForwardEulerTimeStepper) = nothing

#####
##### Checkpointing
#####

# Forward Euler is a self-starting timestepper, so no state needs to be saved
prognostic_state(::ForwardEulerTimeStepper) = nothing
restore_prognostic_state!(ts::ForwardEulerTimeStepper, ::Nothing) = ts
restore_prognostic_state!(ts::ForwardEulerTimeStepper, state) = ts
