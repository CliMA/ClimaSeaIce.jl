using ClimaSeaIce.ThermalBoundaryConditions:
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
            model.ice_thickness,
            Δt,
            grid,
            model.clock,
            model.velocities,
            model.advection,
            model.ice_concentration,
            model.top_surface_temperature,
            model.thermal_boundary_conditions.top,
            model.thermal_boundary_conditions.bottom,
            model.external_thermal_fluxes.top,
            model.internal_thermal_flux,
            model.external_thermal_fluxes.bottom,
            model.ice_consolidation_thickness,
            model.phase_transitions,
            nothing, #model.forcing.h,
            fields(model))

    tick!(model.clock, Δt)

    return nothing
end

@kernel function slab_model_time_step!(ice_thickness, Δt,
                                       grid,
                                       clock,
                                       velocities,
                                       advection,
                                       ice_concentration,
                                       top_temperature,
                                       top_thermal_bc,
                                       bottom_thermal_bc,
                                       top_external_thermal_flux,
                                       internal_thermal_flux,
                                       bottom_external_thermal_flux,
                                       ice_consolidation_thickness,
                                       phase_transitions,
                                       h_forcing,
                                       model_fields)

    i, j = @index(Global, NTuple)

    ℵ = ice_concentration
    h = ice_thickness
    hᶜ = ice_consolidation_thickness
    Qi = internal_thermal_flux
    Qu = top_external_thermal_flux
    Qb = bottom_external_thermal_flux
    Tu = top_temperature
    liquidus = phase_transitions.liquidus

    # Determine top surface temperature
    if !isa(top_thermal_bc, PrescribedTemperature) # update surface temperature?

        consolidated_ice = @inbounds h[i, j, 1] >= hᶜ[i, j, 1]

        if consolidated_ice # slab is consolidated and has an independent surface temperature
            Tu⁻ = @inbounds Tu[i, j, 1]
            Tuⁿ = top_surface_temperature(i, j, grid, top_thermal_bc, Tu⁻, Qi, Qu, clock, model_fields)
        else # slab is unconsolidated and does not have an independent surface temperature
            Tuⁿ = bottom_temperature(i, j, grid, bottom_thermal_bc, liquidus)
        end

        @inbounds Tu[i, j, 1] = Tuⁿ
    end

    Gh = ice_thickness_tendency(i, j, grid, clock,
                                velocities,
                                advection,
                                ice_thickness,
                                ice_consolidation_thickness,
                                top_temperature,
                                bottom_thermal_bc,
                                top_external_thermal_flux,
                                internal_thermal_flux,
                                bottom_external_thermal_flux,
                                phase_transitions,
                                h_forcing,
                                model_fields)

    # Update ice thickness, clipping negative values

    @inbounds begin
        h⁺ = h[i, j, 1] + Δt * Gh
        h[i, j, 1] = max(zero(grid), h⁺)

        # Belongs in update state?
        # That's certainly a simple model for ice concentration
        consolidated_ice = h[i, j, 1] >= hᶜ[i, j, 1]
        ℵ[i, j, 1] = ifelse(consolidated_ice, one(grid), zero(grid))
    end
end

function update_state!(model::SSIM)
    h = model.ice_thickness
    fill_halo_regions!(h)
    return nothing
end

