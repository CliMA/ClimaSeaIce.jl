using Oceananigans.Advection
using ClimaSeaIce.SeaIceMomentumEquations: compute_momentum_tendencies!

function compute_tendencies!(model::SIM, Δt)
    compute_tracer_tendencies!(model)
    compute_momentum_tendencies!(model, model.dynamics, Δt)
    return nothing
end

function compute_tracer_tendencies!(model::SIM)
    grid = model.grid
    arch = architecture(grid)
   
    launch!(arch, grid, :xy,
            _compute_dynamic_ice_variables_tendencies!,
            model.timestepper.Gⁿ,
            grid,
            model.velocities,
            model.advection,
            model.ice_thickness,
            model.ice_concentration)

    # Advance tracers
    for tracer_idx in keys(model.tracers)
        tracer = @inbounds model.tracers[tracer_idx]

        if (tracer_idx ∈ keys(model.timestepper.Gⁿ))
            tendency = @inbounds model.timestepper.Gⁿ[tracer_idx]
            launch!(arch, grid, :xy,
                    _compute_dynamic_tracer_tendency!,
                    tendency,
                    grid,
                    model.velocities,
                    model.advection,
                    tracer)
        end

    end

    return nothing
end

@kernel function _compute_dynamic_tracer_tendency!(Gⁿ, 
                                                   grid,
                                                   velocities,
                                                   advection,
                                                   tracer)
    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) # Assumption! The sea ice is located at the _top_ of the grid

    @inbounds Gⁿ[i, j, 1] = - horizontal_div_Uc(i, j, kᴺ, grid, advection, velocities, tracer)
end

@kernel function _compute_dynamic_ice_variables_tendencies!(Gⁿ, 
                                                            grid,
                                                            velocities,
                                                            advection,
                                                            ice_thickness,
                                                            ice_concentration)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) # Assumption! The sea ice is located at the _top_ of the grid
 
    @inbounds begin
        Gⁿ.h[i, j, 1] = - horizontal_div_Uc(i, j, kᴺ, grid, advection, velocities, ice_thickness)
        Gⁿ.ℵ[i, j, 1] = - horizontal_div_Uc(i, j, kᴺ, grid, advection, velocities, ice_concentration)

        # TODO: BBM rheology needs this!
        # compute_tracer_tendencies!(Gⁿ, i, j, grid, advection, velocities, tracers)
    end
end

const EmptyTuples = Union{NamedTuple{(), Tuple{}}, Tuple{}}

compute_tracer_tendencies!(G, i, j, grid, advection, velocities, ::EmptyTuples) = nothing

function compute_tracer_tendencies!(G, i, j, grid, advection, velocities, tracers)
    # Assumption! The tracer tendencies are the first ones
    for n in keys(tracers)
        if n ∈ keys(G)
            @inbounds G[n][i, j, 1] = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, tracers[n])
        end
    end
end
