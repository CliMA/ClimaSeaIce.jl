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
    water_heat_capacity :: Cp
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

const reference_density = 999.8 # kg m⁻³

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
                                 ice_heat_capacity = 2090.0 / reference_density,
                                 water_heat_capacity = 3991.0 / reference_density,
                                 fusion_enthalpy = 3.3e5 / reference_density,
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
                                 water_heat_capacity,
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
    ϕ = model.state.ϕ
    c = model.ice_heat_capacity
    ℒ = model.fusion_enthalpy

    # set temperature from enthalpy via dH = c dT + ℒ ϕ
    interior(H) .= c .* interior(T) + ℒ * interior(ϕ)

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
    update_diffusivity!(model.closure, model)
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

struct MolecularDiffusivity{C}
    κ_ice :: Float64
    κ_water :: Float64
    κ :: C
end


function MolecularDiffusivity(grid; κ_ice=1e-5, κ_water=1e-6)
    κ = CenterField(grid)
    return MolecularDiffusivity(κ_ice, κ_water, κ)
end

function update_diffusivity!(closure::MolecularDiffusivity, model)
    κ = closure.κ
    κ_ice = closure.κ_ice
    κ_water = closure.κ_water
    ϕ = model.state.ϕ
    grid = model.grid
    arch = grid.architecture
    launch!(arch, grid, :xyz, _compute_molecular_diffusivity!, κ, κ_ice, κ_water, ϕ)
    fill_halo_regions!(κ)
    return nothing
end

@kernel function _compute_molecular_diffusivity!(κ, κ_ice, κ_water, ϕ)
    i, j, k = @index(Global, NTuple) 
    @inbounds begin
        ϕᵢ = ϕ[i, j, k] # porosity or liquid fraction
        κ[i, j, k] = κ_ice * (1 - ϕᵢ) + κ_water * ϕᵢ
    end
end

@inbounds κᶜᶜᶜ_∂zᶜᶜᶠT(i, j, k, grid, κ, T) = ℑzᵃᵃᶠ(i, j, k, grid, κ) * ∂zᶜᶜᶠ(i, j, k, grid, T)
@inline ∂z_κ_∂z_T(i, j, k, grid, closure::MolecularDiffusivity, T) = ∂zᶜᶜᶜ(i, j, k, grid, κᶜᶜᶜ_∂zᶜᶜᶠT, closure.κ, T)

end # module

