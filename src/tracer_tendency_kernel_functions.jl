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

        # for (n, θ) in enumerate(tracers)
        #     @inbounds Gⁿ[n] = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, θ)
        # end
    end
end


