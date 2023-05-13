module ClimaSeaIce

using Oceananigans.BoundaryConditions:
    fill_halo_regions!,
    regularize_field_boundary_conditions,
    FieldBoundaryConditions,
    ValueBoundaryCondition

using Oceananigans.Utils: prettysummary, prettytime
using Oceananigans.Fields: CenterField, ZFaceField, Field, Center, Face, interior
using Oceananigans.Models: AbstractModel
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Utils: launch!
using Oceananigans.Operators

using KernelAbstractions: @kernel, @index

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!

mutable struct ThermodynamicIceModel{Grid,
                                     Tim,
                                     Clk,
                                     Clo,
                                     State,
                                     Cp,
                                     Fu,
                                     Ocean,
                                     Tend} <: AbstractModel{Nothing}
    grid :: Grid
    timestepper :: Tim # unused placeholder for now
    clock :: Clk
    closure :: Clo
    state :: State
    ice_heat_capacity :: Cp
    fusion_enthalpy :: Fu
    ocean_temperature :: Ocean
    tendencies :: Tend
end

const TIM = ThermodynamicIceModel

function Base.show(io::IO, model::TIM)
    clock = model.clock
    print(io, "ThermodynamicIceModel(t=", prettytime(clock.time), ", iteration=", clock.iteration, ")", '\n')
    print(io, "    grid: ", summary(model.grid), '\n')
    print(io, "    ocean_temperature: ", model.ocean_temperature)
end

const reference_density = 1035.0 # kg m⁻³

function default_closure(grid)
    κ = ZFaceField(grid)
    parent(κ) .= 1e-6
    return ScalarIceDiffusivty(κ)
end

"""
    ThermodynamicIceModel(; grid, kw...)

Return a thermodynamic model for ice sandwiched between an atmosphere and ocean.
"""
function ThermodynamicIceModel(; grid,
                                 closure = default_closure(grid),
                                 ice_heat_capacity = 1.0, #2090.0 / reference_density,
                                 fusion_enthalpy = 333.55 / reference_density,
                                 atmosphere_temperature = -10, # ᵒC
                                 ocean_temperature = 0) # ᵒC

    # Build temperature field
    top_T_bc = ValueBoundaryCondition(atmosphere_temperature)
    bottom_T_bc = ValueBoundaryCondition(ocean_temperature)
    T_location = (Center, Center, Center)
    T_bcs = FieldBoundaryConditions(grid, T_location, top=top_T_bc, bottom=bottom_T_bc)
    temperature = CenterField(grid, boundary_conditions=T_bcs)
    enthalpy = CenterField(grid)
    porosity = CenterField(grid)
    state = (T=temperature, H=enthalpy, ϕ=porosity)

    tendencies = (; H=CenterField(grid))
    clock = Clock{eltype(grid)}(0, 0, 1)

    return ThermodynamicIceModel(grid,
                                 nothing,
                                 clock,
                                 closure,
                                 state,
                                 ice_heat_capacity,
                                 fusion_enthalpy,
                                 ocean_temperature,
                                 tendencies)
end

function set!(model::TIM; T=nothing, H=nothing)

    setting_temperature = !isnothing(T)
    setting_enthalpy = !isnothing(H)

    setting_temperature && setting_enthalpy && error("Cannot set both temperature and enthalpy!")

    if setting_temperature
        set!(model.state.T, T)
        update_enthalpy!(model)
    end

    if setting_enthalpy
        set!(model.state.H, H)
        update_temperature!(model)
    end

    return nothing
end

#####
##### Utilities
#####

fields(model::TIM) = model.state
prognostic_fields(model::TIM) = (; model.state.H)

#####
##### Time-stepping
#####

function update_porosity!(model)
    T = model.state.T
    ϕ = model.state.ϕ
    grid = model.grid
    arch = grid.architecture
    launch!(arch, grid, :xyz, _compute_porosity!, ϕ, grid, T)
    return nothing
end

@kernel function _compute_porosity!(ϕ, grid, T)   
    i, j, k = @index(Global, NTuple) 

    FT = eltype(grid)
    Tₘ = zero(FT) # melting temperature

    @inbounds begin
        Tᵢ = T[i, j, k] 
        ϕ[i, j, k] = ifelse(Tᵢ < Tₘ, one(FT), zero(FT))
    end
end

function update_temperature!(model)
    H = model.state.H
    T = model.state.T
    c = model.ice_heat_capacity

    # set temperature from enthalpy via dH = c dT + ℒ ϕ
    interior(T) .= interior(H) ./ c

    arch = model.grid.architecture
    fill_halo_regions!(T, arch, model.clock, fields(model))

    return nothing
end

function update_enthalpy!(model)
    H = model.state.H
    T = model.state.T
    c = model.ice_heat_capacity

    # set temperature from enthalpy via dH = c dT
    interior(H) .= c .* interior(T)

    arch = model.grid.architecture
    fill_halo_regions!(H, arch, model.clock, fields(model))

    return nothing
end

function update_state!(model::TIM)
    grid = model.grid
    arch = grid.architecture
    args = (model.clock, fields(model))
    update_temperature!(model)
    update_porosity!(model)
    return nothing
end

function time_step!(model, Δt; callbacks=nothing)
    grid = model.grid
    arch = grid.architecture
    U = model.state
    G = model.tendencies
    closure = model.closure

    launch!(arch, grid, :xyz, compute_tendencies!, G, grid, closure, U)
    launch!(arch, grid, :xyz, step_fields!, U, G, Δt)
    update_state!(model)
    tick!(model.clock, Δt)

    return nothing
end

@kernel function step_fields!(U, G, Δt)
    i, j, k = @index(Global, NTuple)
    @inbounds U.H[i, j, k] += Δt * G.H[i, j, k]
end

#####
##### Physics
#####

""" Calculate the right-hand-side of the free surface displacement (η) equation. """
@kernel function compute_tendencies!(G, grid, closure, U)
    i, j, k = @index(Global, NTuple)

    # Temperature tendency
    @inbounds G.H[i, j, k] = ∂z_κ_∂z_T(i, j, k, grid, closure, U.T)
end

struct ScalarIceDiffusivity{C}
    κ :: C
end

@inbounds κ_∂zᶜᶜᶠ(i, j, k, grid, κ::Number, T) = κ * ∂zᶜᶜᶠ(i, j, k, grid, T)
@inbounds κ_∂zᶜᶜᶠ(i, j, k, grid, κ::AbstractArray{<:Any, 3}, T) = @inbounds κ[i, j, k] * ∂zᶜᶜᶠ(i, j, k, grid, T)
@inbounds κ_∂zᶜᶜᶠ(i, j, k, grid, κ::AbstractArray{<:Any, 1}, T) = @inbounds κ[k] * ∂zᶜᶜᶠ(i, j, k, grid, T)

@inline ∂z_κ_∂z_T(i, j, k, grid, closure::ScalarIceDiffusivity, T) = ∂zᶜᶜᶜ(i, j, k, grid, κ_∂zᶜᶜᶠ, closure.κ, T)

end # module

