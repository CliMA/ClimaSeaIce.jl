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
            sea_ice_time_step!,
            model.ice_thickness,
            Δt,
            grid,
            model.clock,
            model.velocities,
            model.advection,
            model.ice_concentration,
            model.sea_ice_thermodynamics,
            model.external_heat_fluxes.top,
            model.external_heat_fluxes.bottom,
            nothing, #model.forcing.h,
            fields(model))

    tick!(model.clock, Δt)

    return nothing
end

@kernel function sea_ice_time_step!(ice_thickness, Δt,
                                    grid,
                                    clock,
                                    velocities,
                                    advection,
                                    concentration,
                                    thermodynamics,
                                    top_external_heat_flux,
                                    bottom_external_heat_flux,
                                    h_forcing,
                                    model_fields)

    i, j = @index(Global, NTuple)

    Gh = ice_thickness_tendency(i, j, 1, grid, clock,
                                velocities,
                                advection,
                                ice_thickness,
                                concentration,
                                thermodynamics,
                                top_external_heat_flux,
                                bottom_external_heat_flux,
                                h_forcing,
                                model_fields)
     
    # Update ice thickness, clipping negative values

    @inbounds begin
        h⁺ = h[i, j, 1] + Δt * Gh
        h[i, j, 1] = max(zero(grid), h⁺)
        # Belongs in update state?
        # That's certainly a simple model for ice concentration
        # consolidated_ice = h[i, j, 1] >= hᶜ[i, j, 1]
        # ℵ[i, j, 1] = ifelse(consolidated_ice, one(grid), zero(grid))
    end
end

function update_state!(model::SSIM)
    h = model.ice_thickness
    fill_halo_regions!(h)
    return nothing
end

