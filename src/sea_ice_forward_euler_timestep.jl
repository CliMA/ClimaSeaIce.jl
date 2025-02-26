using ClimaSeaIce.SeaIceMomentumEquations: step_momentum!

const AB2SeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:ForwardEuler}

function time_step!(model::AB2SeaIceModel, Δt; euler=false, callbacks = [])
    
    # Be paranoid and update state at iteration 0
    model.clock.iteration == 0 && update_state!(model)

    compute_tendencies!(model, Δt)
    step_tracers!(model, Δt, 1)

    # TODO: This is an implicit (or split-explicit) step to advance momentum!
    step_momentum!(model, model.dynamics, Δt, 1)

    tick!(model.clock, Δt)
    update_state!(model)

    return nothing
end
