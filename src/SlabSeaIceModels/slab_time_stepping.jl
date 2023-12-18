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

function time_step!(model::SSIM, Δt; callbacks=nothing)
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            slab_model_time_step!,
            model.thickness,
            Δt,
            grid,
            model.clock,
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
            nothing, #model.forcing.h,
            fields(model))

    tick!(model.clock, Δt)

    return nothing
end

@kernel function slab_model_time_step!(thickness, Δt,
                                       grid,
                                       clock,
                                       velocities,
                                       advection,
                                       concentration,
                                       top_temperature,
                                       top_heat_bc,
                                       bottom_heat_bc,
                                       top_external_heat_flux,
                                       internal_heat_flux,
                                       bottom_external_heat_flux,
                                       consolidation_thickness,
                                       phase_transitions,
                                       forcing,
                                       model_fields)

    i, j = @index(Global, NTuple)

    ℵ  = concentration
    h  = thickness
    hᶜ = consolidation_thickness
    Qi = internal_heat_flux
    Qu = top_external_heat_flux
    Qb = bottom_external_heat_flux
    Tu = top_temperature
    liquidus = phase_transitions.liquidus

    # Determine top surface temperature
    if !isa(top_heat_bc, PrescribedTemperature) # update surface temperature?

        consolidated_ice = @inbounds h[i, j, 1] >= hᶜ[i, j, 1]

        if consolidated_ice # slab is consolidated and has an independent surface temperature
            Tu⁻ = @inbounds Tu[i, j, 1]
            Tuⁿ = top_surface_temperature(i, j, grid, top_heat_bc, Tu⁻, Qi, Qu, clock, model_fields)
        else # slab is unconsolidated and does not have an independent surface temperature
            Tuⁿ = bottom_temperature(i, j, grid, bottom_heat_bc, liquidus)
        end

        @inbounds Tu[i, j, 1] = Tuⁿ
    end

    Gh = thickness_tendency(i, j, grid, clock,
                            velocities,
                            advection,
                            thickness,
                            concentration,
                            consolidation_thickness,
                            top_temperature,
                            bottom_heat_bc,
                            top_external_heat_flux,
                            internal_heat_flux,
                            bottom_external_heat_flux,
                            phase_transitions,
                            forcing,
                            model_fields)

    Gℵ = concentration_tendency(i, j, grid, clock,
                                velocities,
                                advection,
                                thickness,
                                concentration,
                                consolidation_thickness,
                                top_temperature,
                                bottom_heat_bc,
                                top_external_heat_flux,
                                internal_heat_flux,
                                bottom_external_heat_flux,
                                phase_transitions,
                                forcing,
                                model_fields)

    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁺ = h[i, j, 1] + Δt * Gh
        h[i, j, 1] = max(zero(grid), h⁺)

        # Belongs in update state?
        # That's certainly a simple model for ice concentration
        consolidated_ice = h[i, j, 1] >= hᶜ[i, j, 1]
        ℵ⁺ = ℵ[i, j, 1] + Δt * Gℵ
        ℵ[i, j, 1] = ifelse(consolidated_ice, max(zero(grid), ℵ⁺), zero(grid))
    end
    
end

function update_state!(model::SSIM)
    h = model.thickness
    ℵ = model.concentration
    fill_halo_regions!(h)
    fill_halo_regions!(ℵ)
    return nothing
end

