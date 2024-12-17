const RK3SeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:RungeKutta3TimeStepper}

function time_step!(model::RK3SeaIceModel, Δt)

    # Be paranoid and update state at iteration 0, in case run! is not used:
    model.clock.iteration == 0 && update_state!(model, callbacks)

    γ¹ = model.timestepper.γ¹
    γ² = model.timestepper.γ²
    γ³ = model.timestepper.γ³

    ζ² = model.timestepper.ζ²
    ζ³ = model.timestepper.ζ³

    first_stage_Δt  = γ¹ * Δt
    second_stage_Δt = (γ² + ζ²) * Δt
    third_stage_Δt  = (γ³ + ζ³) * Δt

    #
    # First stage
    #

    compute_tracer_tendencies!(model)
    step_tracers!(model, Δt, 1)

    # TODO: This is an implicit (or split-explicit) step to advance momentum!
    # do we need to pass Δt here or first_stage_Δt?
    step_momentum!(model, model.ice_dynamics, Δt)
    store_tendencies!(model)

    tick!(model.clock, first_stage_Δt)
    update_state!(model)

    #
    # Second stage
    #

    compute_tracer_tendencies!(model)
    step_tracers!(model, Δt, 2)

    # TODO: This is an implicit (or split-explicit) step to advance momentum!
    # do we need to pass Δt here or second_stage_Δt?
    step_momentum!(model, model.ice_dynamics, Δt)
    store_tendencies!(model)

    tick!(model.clock, second_stage_Δt)
    update_state!(model)

    #
    # Third stage
    #

    compute_tracer_tendencies!(model)
    step_tracers!(model, Δt, 3)

    # TODO: This is an implicit (or split-explicit) step to advance momentum!
    # do we need to pass Δt here or third_stage_Δt?
    step_momentum!(model, model.ice_dynamics, Δt)
    store_tendencies!(model)

    tick!(model.clock, third_stage_Δt)
    update_state!(model)

    return nothing
end

function timestepping_coefficients(ts::RungeKutta3TimeStepper, substep) 
    if substep == 1
        return ts.γ¹, zero(ts.γ¹)
    elseif substep == 2
        return ts.γ², ts.ζ²
    elseif substep == 3
        return ts.γ³, ts.ζ³
    end
end