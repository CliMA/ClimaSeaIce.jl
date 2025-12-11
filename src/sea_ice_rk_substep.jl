using Oceananigans.TimeSteppers: SplitRungeKuttaTimeStepper
using ClimaSeaIce.SeaIceDynamics: time_step_momentum!
using ClimaSeaIce.SeaIceThermodynamics: thermodynamic_time_step!

import Oceananigans.TimeSteppers: rk_substep!, cache_previous_fields!, step_lagrangian_particles!

const RKSeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:SplitRungeKuttaTimeStepper}

step_lagrangian_particles!(model::SeaIceModel, Δτ) = nothing

function cache_previous_fields!(model::RKSeaIceModel)
    previous_fields = model.timestepper.Ψ⁻
    model_fields = prognostic_fields(model)
    grid = model.grid
    arch = architecture(grid)

    for name in keys(model_fields)
        Ψ⁻ = previous_fields[name]
        Ψⁿ = model_fields[name]
        parent(Ψ⁻) .= parent(Ψⁿ)
    end

    return nothing
end

# We separate the thermodynamic step from the advection (dynamic) step.
# The thermodynamic step is column physics and is performed all at once.
function rk_substep!(model::RKSeaIceModel, Δτ, callbacks)

    thermodynamic_time_step!(model, model.ice_thermodynamics, Δτ)

    # Compute advective tendencies and update advected tracers
    compute_tendencies!(model, Δτ)
    dynamic_time_step!(model, Δτ)

    # This is an implicit (or split-explicit) step to advance momentum.
    time_step_momentum!(model, model.dynamics, Δτ)

    return nothing
end
