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
import Oceananigans.TimeSteppers: compute_tendencies!, time_step!, ab2_step!, store_tendencies!

function compute_tendencies!(model::SSIM)
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_slab_model_tendencies!,
            model.timestepper.Gⁿ,
            model.thickness,
            Δt,
            grid,
            model.clock,
            model.velocities,
            model.ocean_velocities,
            model.advection,
            model.coriolis,
            model.concentration,
            model.top_surface_temperature,
            model.heat_boundary_conditions.top,
            model.heat_boundary_conditions.bottom,
            model.external_momentum_stress.u.top,
            model.external_momentum_stress.u.bottom,
            model.external_momentum_stress.v.top,
            model.external_momentum_stress.v.bottom,
            model.external_heat_fluxes.top,
            model.internal_heat_flux,
            model.external_heat_fluxes.bottom,
            model.consolidation_thickness,
            model.phase_transitions,
            nothing, #model.forcing
            fields(model))

    return nothing
end

function ab2_step!(model::SSIM, Δt, χ)
    grid = model.grid
    arch = architecture(grid)

    h  = model.thickness
    hᶜ = model.consolidation_thickness
    ℵ  = model.concentration
    u, v = model.velocities

    launch!(arch, grid, :xyz, _ab2_step_slab_model!, h, ℵ, u, v, model.timestepper.Gⁿ, model.timestepper.G⁻, hᶜ, Δt, χ)

    return nothing
end

store_tendencies!(model::SSIM) = 
    launch!(architecture(model.grid), model.grid, :xyz, _store_all_tendencies!, 
                         model.timestepper.G⁻, model.timestepper.Gⁿ, Val(length(model.timestepper.Gⁿ)))

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

    compute_tendencies!(model; callbacks)
    ab2_step!(model, Δt, χ)
    store_tendencies!(model)

    tick!(model.clock, Δt)
    update_state(model)

    return nothing
end

function update_state!(model::SSIM)
    fields = prognostic_fields(model)
    fill_halo_regions!(fields)
    return nothing
end

@kernel function _ab2_step_slab_model!(h, ℵ, u, v, Gⁿ, G⁻, hᶜ, Δt, χ)
    i, j = @index(Global, NTuple)

    FT = eltype(χ)
    Δt = convert(FT, Δt)
    one_point_five = convert(FT, 1.5)
    oh_point_five  = convert(FT, 0.5)

    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁺ = h[i, j, 1] + Δt * ((one_point_five + χ) * Gⁿ.h[i, j, 1] - (oh_point_five + χ) * G⁻.h[i, j, 1])
        h[i, j, 1] = max(zero(grid), h⁺)

        # Belongs in update state?
        # That's certainly a simple model for ice concentration
        consolidated_ice = h[i, j, 1] >= hᶜ[i, j, 1]
        ℵ⁺ = ℵ[i, j, 1] + Δt * ((one_point_five + χ) * Gⁿ.ℵ[i, j, 1] - (oh_point_five + χ) * G⁻.ℵ[i, j, 1])
        ℵ[i, j, 1] = ifelse(consolidated_ice, max(zero(grid), ℵ⁺), zero(grid))

        u[i, j, 1] += Δt * ((one_point_five + χ) * Gⁿ.u[i, j, 1] - (oh_point_five + χ) * G⁻.u[i, j, 1])
        v[i, j, 1] += Δt * ((one_point_five + χ) * Gⁿ.v[i, j, 1] - (oh_point_five + χ) * G⁻.v[i, j, 1])
    end 
end

@kernel function _store_all_tendencies!(G⁻, Gⁿ, ::Val{N}) where N
    i, j = @index(Global, NTuple)

    @unroll for n in 1:N
        G⁻[n][i, j, 1] = Gⁿ[n][i, j, 1]
    end
end