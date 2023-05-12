module ClimaSeaIce

using Oceananigans.BoundaryConditions:
    fill_halo_regions!,
    regularize_field_boundary_conditions,
    FieldBoundaryConditions
    ValueBoundaryCondition

using Oceananigans.Fields: CenterField, Field, Center, Face
using Oceananigans.Models: AbstractModel
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Utils: launch!
using Oceananigans.Operators

using KernelAbstractions: @kernel, @index

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!

mutable struct ThermodynamicIceModel{Grid, Tim, Clk, Clo, Temp, Ocean, Tend} <: AbstractModel{Nothing}
    grid :: Grid
    timestepper :: Tim # unused placeholder for now
    clock :: Clk
    closure :: Clo
    temperature :: Temp
    ocean_temperature :: Ocean
    tendencies :: Tend
end

const TIM = ThermodynamicIceModel

"""
    ThermodynamicIceModel(; grid, kw...)

Return a thermodynamic model for ice sandwiched between an atmosphere and ocean.
"""
function ThermodynamicIceModel(; grid,
                                 closure = ConductiveIceModel(1e-6),
                                 atmosphere_temperature = -10, # ᵒC
                                 ocean_temperature = 0) # ᵒC

    # Build temperature field
    top_T_bc = ValueBoundaryCondition(atmosphere_temperature)
    T_location = (Center, Center, Center)
    T_bcs = FieldBoundaryConditions(grid, T_location, top=top_T_bc)
    temperature = CenterField(grid, boundary_conditions=T_bcs)

    tendencies = (; T=CenterField(grid))
    clock = Clock{eltype(grid)}(0, 0, 1)

    return ThermodynamicIceModel(grid,
                                   nothing,
                                   clock,
                                   closure,
                                   temperature,
                                   ocean_temperature,
                                   tendencies)
end

#####
##### Utilities
#####

fields(model::TIM) = (; T=model.temperature)
prognostic_fields(model::TIM) = fields(model)

#####
##### Time-stepping
#####

function update_state!(model::TIM)
    grid = model.grid
    arch = grid.architecture
    args = (model.clock, fields(model))
    fill_halo_regions!(model.temperature, arch, model.clock, fields(model))
    return nothing
end

function time_step!(model, Δt; callbacks=nothing)
    grid = model.grid
    arch = grid.architecture
    T = model.temperature
    ∂t_T = model.tendencies.T
    closure = model.closure

    launch!(arch, grid, :xyz, calculate_tendencies!, ∂t_T, grid, closure, T)
    launch!(arch, grid, :xyz, step_fields!, T, ∂t_T, Δt)
    update_state!(model)
    tick!(model.clock, Δt)

    return nothing
end

@kernel function step_fields!(T, ∂t_T, Δt)
    i, j, k = @index(Global, NTuple)
    @inbounds T[i, j, k] += Δt * ∂t_T[i, j, k]
end

#####
##### Physics
#####

""" Calculate the right-hand-side of the free surface displacement (η) equation. """
@kernel function calculate_tendencies!(∂t_T, grid, closure, T)
    i, j, k = @index(Global, NTuple)

    # Temperature tendency
    @inbounds ∂t_T[i, j, k] = ∂z_κ_∂z_T(i, j, k, grid, closure, T)
end

struct ConstantIceDiffusivity{C}
    κ :: C
end

@inline ∂z_κ_∂z_T(i, j, k, grid, closure::ConstantIceDiffusivity, T) =
    closure.κ * ∂zᶜᶜᶜ(i, j, k, grid, ∂zᶜᶜᶠ, T)

end # module

