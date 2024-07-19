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

mask_immersed_field!(::ConstantField) = nothing
mask_immersed_field!(::ZeroField)     = nothing

function time_step!(model::SIM, Δt; callbacks=nothing, euler=false)
    
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
    ab2_timestepper.χ = χ

    # Ensure zeroing out all previous tendency fields to avoid errors in
    # case G⁻ includes NaNs. See https://github.com/CliMA/Oceananigans.jl/issues/2259
    if euler
        @debug "Taking a forward Euler step."
        for field in ab2_timestepper.G⁻
            !isnothing(field) && fill!(field, 0)
        end
    end

    # TODO: This is an implicit (or split-explicit) step to advance momentum!
    step_momentum!(model, model.ice_dynamics, Δt, χ)

    compute_tracer_tendencies!(model; callbacks)
    ab2_step_tracers!(model, Δt, χ)

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

function ab2_step_tracers!(model::SIM, Δt, χ)
    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ℵ  = model.ice_concentration
    tracers = model.tracers

    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻

    launch!(arch, grid, :xyz, _ab2_step_tracers!, h, ℵ, tracers, Gⁿ, G⁻, Δt, χ)

    return nothing
end

function store_tendencies!(model::SIM) 

    grid = model.grid
    arch = architecture(grid)
    Nx, Ny, _ = size(grid)

    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻
    Nt = length(Gⁿ)

    params = KernelParameters((Nx, Ny, Nt), (0, 0, 0))
    launch!(arch, model.grid, params, _store_all_tendencies!, G⁻, Gⁿ)

    return nothing
end

# Thickness and concentration are updated using an AB2 scheme
# We compute hⁿ⁺¹ and ℵⁿ⁺¹ in the same kernel to account for ridging: 
# if ℵ > 1, we reset the concentration to 1 and adjust the thickness 
# to conserve the total ice volume in the cell.
@kernel function _ab2_step_tracers!(h, ℵ, tracers, Gⁿ, G⁻, Δt, χ)
    i, j, k = @index(Global, NTuple)

    FT = eltype(χ)
    Δt = convert(FT, Δt)
    one_point_five = convert(FT, 1.5)
    oh_point_five  = convert(FT, 0.5)

    Ghⁿ = Gⁿ.h
    Gℵⁿ = Gⁿ.ℵ
    
    Gh⁻ = G⁻.h
    Gℵ⁻ = G⁻.ℵ

    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁺ = h[i, j, k] + Δt * ((one_point_five + χ) * Ghⁿ[i, j, k] - (oh_point_five + χ) * Gh⁻[i, j, k])
        h⁺ = max(zero(FT), h⁺)

        # Belongs in update state?
        # That's certainly a simple model for ice concentration
        ℵ⁺ = ℵ[i, j, k] + Δt * ((one_point_five + χ) * Gℵⁿ[i, j, k] - (oh_point_five + χ) * Gℵ⁻[i, j, k])
        ℵ⁺ = max(zero(FT), ℵ⁺)
        
        # Ridging! if ℵ > 1, we reset the concentration to 1 and increase the thickness accordingly
        # to maintain a constant ice volume
        h[i, j, k] = ifelse(ℵ⁺ > 1, h⁺ * ℵ⁺, h⁺)
        ℵ[i, j, k] = ifelse(ℵ⁺ > 1, one(FT), ℵ⁺)
    end 
end

@kernel function _store_all_tendencies!(G⁻, Gⁿ) 
    i, j, n = @index(Global, NTuple)
    @inbounds G⁻[n][i, j, 1] = Gⁿ[n][i, j, 1]
end

function update_state!(model::SIM, callbacks = nothing)
    
    foreach(prognostic_fields(model)) do field
        mask_immersed_field!(field)
    end

    fill_halo_regions!(prognostic_fields(model), model.clock, fields(model))

    return nothing
end
