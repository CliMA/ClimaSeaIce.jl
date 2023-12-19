using ClimaSeaIce.HeatBoundaryConditions:
    PrescribedTemperature,
    bottom_temperature,
    top_surface_temperature,
    getflux

using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.TimeSteppers: tick!
using Oceananigans.Utils: launch!

using KernelAbstractions: @index, @kernel
using Oceananigans.TimeSteppers: ab2_step_field!
import Oceananigans.TimeSteppers: time_step!

function compute_tracer_tendencies!(model::SSIM; callbacks = nothing)
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_tracer_tendencies!,
            model.timestepper.Gⁿ,
            grid,
            model.clock,
            model.thickness,
            model.concentration,
            model.velocities,
            model.advection,
            model.concentration,
            model.top_surface_temperature,
            model.heat_boundary_conditions.top,
            model.heat_boundary_conditions.bottom,
            model.external_heat_fluxes.top,
            model.internal_heat_flux,
            model.external_heat_fluxes.bottom,
            model.consolidation_thickness,
            model.phase_transitions,
            nothing, #model.forcing
            fields(model))

    return nothing
end

function ab2_step_tracers!(model::SSIM, Δt, χ)
    grid = model.grid
    arch = architecture(grid)

    h  = model.thickness
    hᶜ = model.consolidation_thickness
    ℵ  = model.concentration

    Ghⁿ = model.timestepper.Gⁿ.h
    Gh⁻ = model.timestepper.G⁻.h
    Gℵⁿ = model.timestepper.Gⁿ.ℵ
    Gℵ⁻ = model.timestepper.G⁻.ℵ    

    launch!(arch, grid, :xyz, _ab2_step_tracers!, h, ℵ, Ghⁿ, Gℵⁿ, Gh⁻, Gℵ⁻, hᶜ, Δt, χ)

    return nothing
end

store_tracer_tendencies!(model::SSIM) = 
    launch!(architecture(model.grid), model.grid, :xyz, _store_all_tendencies!, 
                         model.timestepper.G⁻[2:end], model.timestepper.Gⁿ[2:end], Val(length(model.timestepper.Gⁿ[2:end])))

function time_step!(model::SSIM, Δt; callbacks=nothing, euler=false)

    χ = ifelse(euler, convert(eltype(model.grid), -0.5), model.timestepper.χ)

    if euler
        @debug "Taking a forward Euler step."
        # Ensure zeroing out all previous tendency fields to avoid errors in
        # case G⁻ includes NaNs. See https://github.com/CliMA/Oceananigans.jl/issues/2259
        for field in model.timestepper.G⁻
            !isnothing(field) && fill!(field, 0)
        end
    end

    compute_tracer_tendencies!(model; callbacks)
    ab2_step_tracers!(model, Δt, χ)
    store_tracer_tendencies!(model)

    # TODO: Add the rheology here!
    advance_momentum!(model, Δt, χ)

    tick!(model.clock, Δt)
    update_state(model)

    return nothing
end

function update_state!(model::SSIM)
    fields = prognostic_fields(model)
    fill_halo_regions!(fields)
    return nothing
end

@kernel function _ab2_step_slab_model!(h, ℵ, Ghⁿ, Gℵⁿ, Gh⁻, Gℵ⁻, hᶜ, Δt, χ)
    i, j = @index(Global, NTuple)

    FT = eltype(χ)
    Δt = convert(FT, Δt)
    one_point_five = convert(FT, 1.5)
    oh_point_five  = convert(FT, 0.5)

    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁺ = h[i, j, 1] + Δt * ((one_point_five + χ) * Ghⁿ[i, j, 1] - (oh_point_five + χ) * Gh⁻[i, j, 1])
        h[i, j, 1] = max(zero(grid), h⁺)

        # Belongs in update state?
        # That's certainly a simple model for ice concentration
        consolidated_ice = h[i, j, 1] >= hᶜ[i, j, 1]
        ℵ⁺ = ℵ[i, j, 1] + Δt * ((one_point_five + χ) * Gℵⁿ[i, j, 1] - (oh_point_five + χ) * Gℵ⁻[i, j, 1])
        ℵ[i, j, 1] = ifelse(consolidated_ice, max(zero(grid), ℵ⁺), zero(grid))
    end 
end

@kernel function _store_all_tendencies!(G⁻, Gⁿ, ::Val{N}) where N
    i, j = @index(Global, NTuple)

    @unroll for n in 1:N
        G⁻[n][i, j, 1] = Gⁿ[n][i, j, 1]
    end
end

# This will be the implicit step whenever rheology is implemented!
function advance_momentum!(model, Δt, χ)
    grid = model.grid
    arch = architecture(grid)

    u, v = model.velocities
    Guⁿ = model.timestepper.Gⁿ.u
    Gu⁻ = model.timestepper.G⁻.u
    Gvⁿ = model.timestepper.Gⁿ.v
    Gv⁻ = model.timestepper.G⁻.v

    launch!(arch, grid, :xyz, 
            _compute_momentum_tendencies!, 
            model.timestepper.Gⁿ,
            grid,
            model.clock,
            model.velocities,
            model.ocean_velocities,
            model.coriolis,
            model.external_momentum_stress.u.top,
            model.external_momentum_stress.u.bottom,
            model.external_momentum_stress.v.top,
            model.external_momentum_stress.v.bottom,
            nothing, #model.forcing
            fields(model))

    launch!(arch, grid, :xyz, ab2_step_field!, u, Δt, χ, Guⁿ, Gu⁻)
    launch!(arch, grid, :xyz, ab2_step_field!, v, Δt, χ, Gvⁿ, Gv⁻)

    return nothing
end