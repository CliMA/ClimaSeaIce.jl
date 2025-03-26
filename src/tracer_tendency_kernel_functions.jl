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
            _compute_dynamic_tracer_tendencies!,
            model.timestepper.Gⁿ,
            grid,
            model.velocities,
            model.advection,
            model.ice_thickness,
            model.ice_concentration,
            model.tracers)

    return nothing
end

@kernel function _compute_dynamic_tracer_tendencies!(Gⁿ, 
                                                     grid,
                                                     velocities,
                                                     advection,
                                                     ice_thickness,
                                                     ice_concentration,
                                                     tracers)

    i, j = @index(Global, NTuple)
    
    @inbounds begin
        Gⁿ.h[i, j, 1] = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, ice_thickness)
        Gⁿ.ℵ[i, j, 1] = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, ice_concentration)

        # TODO: BBM rheology needs this!
        # compute_tracer_tendencies!(Gⁿ, i, j, grid, advection, velocities, tracers)
    end
end

# const EmptyTuples = Union{NamedTuple{(), Tuple{}}, Tuple{}}

# compute_tracer_tendencies!(G, i, j, grid, advection, velocities, ::EmptyTuples) = nothing

# function compute_tracer_tendencies!(G, i, j, grid, advection, velocities, tracers)
#     # Assumption! The tracer tendencies are the first ones
#     for n in eachindex(G)
#         @inbounds G[n][i, j, 1] = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, tracers[n])
#     end
# end
