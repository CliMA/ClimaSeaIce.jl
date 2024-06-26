using ClimaSeaIce.HeatBoundaryConditions:
    PrescribedTemperature,
    bottom_temperature,
    top_surface_temperature,
    getflux

using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.TimeSteppers: tick!
using Oceananigans.Utils: launch!, KernelParameters
import Oceananigans.ImmersedBoundaries: mask_immersed_field!

using Oceananigans: prognostic_fields
using KernelAbstractions: @index, @kernel
using Oceananigans.TimeSteppers: ab2_step_field!
import Oceananigans.TimeSteppers: time_step!

mask_immersed_field!(field::ConstantField) = nothing
mask_immersed_field!(field::ZeroField)     = nothing

function compute_tracer_tendencies!(model::SSIM; callbacks = nothing)
    grid = model.grid
    arch = architecture(grid)
   
    launch!(arch, grid, :xyz,
            _compute_tracer_tendencies!,
            model.timestepper.Gⁿ,
            grid,
            model.clock,
            model.ice_thickness,
            model.concentration,
            model.velocities,
            model.advection,
            model.top_surface_temperature,
            model.heat_boundary_conditions.top,
            model.heat_boundary_conditions.bottom,
            model.external_heat_fluxes.top,
            model.internal_heat_flux,
            model.external_heat_fluxes.bottom,
            model.ice_consolidation_thickness,
            model.phase_transitions,
            nothing, #model.forcing
            fields(model))

    return nothing
end

function ab2_step_tracers!(model::SSIM, Δt, χ)
    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    hᶜ = model.ice_consolidation_thickness
    ℵ  = model.concentration

    Ghⁿ = model.timestepper.Gⁿ.h
    Gh⁻ = model.timestepper.G⁻.h
    Gℵⁿ = model.timestepper.Gⁿ.ℵ
    Gℵ⁻ = model.timestepper.G⁻.ℵ    

    launch!(arch, grid, :xyz, _ab2_step_tracers!, h, ℵ, Ghⁿ, Gℵⁿ, Gh⁻, Gℵ⁻, hᶜ, Δt, χ)

    return nothing
end

function store_tendencies!(model::SSIM) 

    grid = model.grid
    arch = architecture(grid)
    Nx, Ny, _ = size(grid)

    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻
    Nt = length(Gⁿ)

    params = KernelParameters((Nx, Ny, Nt), (0, 0, 0))
    launch!(architecture(model.grid), model.grid, params, _store_all_tendencies!, G⁻, Gⁿ)

    return nothing
end

function update_state!(model::SSIM)
    
    foreach(prognostic_fields(model)) do field
        mask_immersed_field!(field)
    end

    fill_halo_regions!(prognostic_fields(model), model.clock, fields(model))

    return nothing
end

function time_step!(model::SSIM, Δt; callbacks=nothing, euler=false)

    # Shenanigans for properly starting the AB2 loop with an Euler step
    euler = euler || (Δt != model.timestepper.previous_Δt)
    
    χ = ifelse(euler, convert(eltype(model.grid), -0.5), model.timestepper.χ)

    if euler
        @debug "Taking a forward Euler step."
        # Ensure zeroing out all previous tendency fields to avoid errors in
        # case G⁻ includes NaNs. See https://github.com/CliMA/Oceananigans.jl/issues/2259
        for field in model.timestepper.G⁻
            !isnothing(field) && fill!(field, 0)
        end
    end

    model.timestepper.previous_Δt = Δt
    
    compute_tracer_tendencies!(model; callbacks)
    ab2_step_tracers!(model, Δt, χ)

    # TODO: This is an implicit (or split-explicit) step to advance momentum!
    step_momentum!(model, model.momentum_solver, Δt, χ)

    # Only the tracers are advanced through an AB2 scheme,
    # so only tracers' tendencies are stored
    store_tendencies!(model)

    tick!(model.clock, Δt)
    update_state!(model)

    return nothing
end

# Thickness and concentration are updated using an AB2 scheme
# We compute hⁿ⁺¹ and ℵⁿ⁺¹ in the same kernel to account for ridging: 
# if ℵ > 1, we reset the concentration to 1 and adjust the thickness 
# to conserve the total ice volume in the cell.
@kernel function _ab2_step_tracers!(h, ℵ, Ghⁿ, Gℵⁿ, Gh⁻, Gℵ⁻, hᶜ, Δt, χ)
    i, j = @index(Global, NTuple)

    FT = eltype(χ)
    Δt = convert(FT, Δt)
    one_point_five = convert(FT, 1.5)
    oh_point_five  = convert(FT, 0.5)

    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁺ = h[i, j, 1] + Δt * ((one_point_five + χ) * Ghⁿ[i, j, 1] - (oh_point_five + χ) * Gh⁻[i, j, 1])
        h[i, j, 1] = max(0, h⁺)

        # Belongs in update state?
        # That's certainly a simple model for ice concentration
        consolidated_ice = h[i, j, 1] >= hᶜ[i, j, 1]
        ℵ⁺ = ℵ[i, j, 1] + Δt * ((one_point_five + χ) * Gℵⁿ[i, j, 1] - (oh_point_five + χ) * Gℵ⁻[i, j, 1])
        ℵ[i, j, 1] = ifelse(consolidated_ice, max(0, ℵ⁺), 0)
        
        # Ridging! if ℵ > 1, we reset the concentration to 1 and increase the thickness accordingly
        # to maintain a constant ice volume
        h[i, j, 1] = ifelse(ℵ[i, j, 1] > 1, h[i, j, 1] * ℵ[i, j, 1], h[i, j, 1])
        ℵ[i, j, 1] = ifelse(ℵ[i, j, 1] > 1, 1, ℵ[i, j, 1])
    end 
end

@kernel function _store_all_tendencies!(G⁻, Gⁿ) where N
    i, j, n = @index(Global, NTuple)
    @inbounds G⁻[n][i, j, 1] = Gⁿ[n][i, j, 1]
end

# Fallback for no sea - ice dynamics
step_momentum!(args...) = nothing

