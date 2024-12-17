using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.TimeSteppers: tick!
using Oceananigans.Utils: launch!, KernelParameters
import Oceananigans.ImmersedBoundaries: mask_immersed_field!

using Oceananigans: prognostic_fields
using KernelAbstractions: @index, @kernel
using Oceananigans.TimeSteppers: ab2_step_field!
import Oceananigans.TimeSteppers: time_step!

using ClimaSeaIce.SeaIceDynamics: step_momentum!

const AB2SeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:QuasiAdamsBashforth2TimeStepper}

function time_step!(model::AB2SeaIceModel, Δt; callbacks=nothing, euler=false)
    
    # Be paranoid and update state at iteration 0
    model.clock.iteration == 0 && update_state!(model, callbacks)

    ab2_timestepper = model.timestepper

    # Change the default χ if necessary, which occurs if:
    #   * We detect that the time-step size has changed.
    #   * We detect that this is the "first" time-step, which means we
    #     need to take an euler step. Note that model.clock.last_Δt is
    #     initialized as Inf
    #   * The user has passed euler=true to time_step!
    euler = euler || (Δt != model.clock.last_Δt)
    
    # If euler, then set χ = -0.5
    minus_point_five = convert(eltype(model.grid), -0.5)
    χ = ifelse(euler, minus_point_five, ab2_timestepper.χ)

    # Set time-stepper χ (this is used in ab2_step!, but may also be used elsewhere)
    χ₀ = ab2_timestepper.χ # Save initial value
    ab2_timestepper.χ = χ

    # Ensure zeroing out all previous tendency fields to avoid errors in
    # case G⁻ includes NaNs. See https://github.com/CliMA/Oceananigans.jl/issues/2259
    if euler
        @debug "Taking a forward Euler step."
        for field in ab2_timestepper.G⁻
            !isnothing(field) && fill!(field, 0)
        end
    end

    compute_tracer_tendencies!(model; callbacks)
    step_tracers!(model, Δt, 1)

    # TODO: This is an implicit (or split-explicit) step to advance momentum!
    step_momentum!(model, model.ice_dynamics, Δt)

    # Only the tracers are advanced through an AB2 scheme 
    # (velocities are stepped in the dynamics step)
    # so only tracers' tendencies are stored
    store_tendencies!(model)

    tick!(model.clock, Δt)
    update_state!(model)

    return nothing
end

function compute_tracer_tendencies!(model::SIM; callbacks = nothing)
    grid = model.grid
    arch = architecture(grid)
   
    launch!(arch, grid, :xyz,
            _compute_tracer_tendencies!,
            model.timestepper.Gⁿ,
            model.ice_thickness,
            grid,
            model.clock,
            model.velocities,
            model.advection,
            model.ice_concentration,
            model.ice_thermodynamics,
            model.external_heat_fluxes.top,
            model.external_heat_fluxes.bottom,
            nothing, #model.forcing.h,
            fields(model))

    return nothing
end

function timestepping_coefficients(ts::QuasiAdamsBashforth2TimeStepper, args...) 
    χ  = ts.χ
    FT = eltype(χ)
    α  = + convert(FT, 1.5) + χ
    β  = - convert(FT, 0.5) + χ
    return α, β
end
