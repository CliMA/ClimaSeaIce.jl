using Oceananigans.Advection
using ClimaSeaIce.SeaIceDynamics: compute_momentum_tendencies!

function compute_tendencies!(model::SIM, Δt)
    compute_tracer_tendencies!(model)
    compute_momentum_tendencies!(model, model.dynamics, Δt)
    return nothing
end

function compute_tracer_tendencies!(model::SIM)
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy,
            _compute_dynamic_tracer_tendencies!,
            model.timestepper.Gⁿ,
            grid,
            model.velocities,
            model.advection,
            model.ice_thickness,
            model.ice_concentration,
            model.snow_thickness,
            model.tracers)

    return nothing
end

@kernel function _compute_dynamic_tracer_tendencies!(Gⁿ,
                                                     grid,
                                                     velocities,
                                                     advection,
                                                     ice_thickness,
                                                     ice_concentration,
                                                     snow_thickness,
                                                     tracers)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) # Assumption! The sea ice is located at the _top_ of the grid

    @inbounds begin
        Gⁿ.h[i, j, 1] = - horizontal_div_Uc(i, j, kᴺ, grid, advection, velocities, ice_thickness)
        Gⁿ.ℵ[i, j, 1] = - horizontal_div_Uc(i, j, kᴺ, grid, advection, velocities, ice_concentration)
    end

    compute_snow_advection_tendency!(i, j, kᴺ, Gⁿ, grid, advection, velocities, ice_concentration, snow_thickness)
end

@inline compute_snow_advection_tendency!(i, j, k, Gⁿ, grid, advection, velocities, ℵ, ::Nothing) = nothing

@inline function compute_snow_advection_tendency!(i, j, k, Gⁿ, grid, advection, velocities, ℵ, hs)
    @inbounds Gⁿ.hs[i, j, 1] = - div_Uℵh(i, j, k, grid, advection, velocities, ℵ, hs)
    return nothing
end
