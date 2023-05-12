module ClimaSeaIce

using Oceananigans
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization

# Simulations interface
import Oceananigans: fields
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!

struct DiffusiveIceModel <: AbstractModel
    grid
    clock
    closure
    temperature
    tendencies
end

function DiffusiveIceModel(; grid, diffusivity)
    vitd = VerticallyImplicitTimeDiscretization()
    closure = VerticalScalarDiffusivity(vitd, κ=κ)
    temperature = CenterField(grid)
    tendencies = (; T=CenterField(grid))
    return DiffusiveIceModel(grid, closure, temperature, tendencies)
end

@kernel function step_fields!(T, ∂T∂t, Δt)
    i, j, k = @index(Global, NTuple)
    @inbounds begin
        T[i, j, k] += Δt * ∂t_T[i, j, k]
    end
end

function time_step!(model, Δt; callbacks=nothing)
    calculate_tendencies!(model)

    arch = model.grid.architecture
    launch!(arch, model.grid, :xyz,
            step_fields!, model.temperature, model.tendencies.T, Δt)

    update_state!(model)
    tick!(model.clock, Δt)

    return nothing
end

""" Calculate the right-hand-side of the free surface displacement (η) equation. """
@kernel function _calculate_interior_tendencies!(∂t_T, grid, κ, T)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        ∂t_T = κ * ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, T)
    end
end

function calculate_tendencies!(∂t_T, grid, κ, T)
    arch = grid.architecture

    launch!(arch, grid, :xyz,
            _calculate_interior_tendencies!, ∂t_T, grid, κ, T)

    return nothing
end

end # module
