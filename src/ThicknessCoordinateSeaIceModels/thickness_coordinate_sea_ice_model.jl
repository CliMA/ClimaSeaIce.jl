using Oceananigans.BoundaryConditions:
    fill_halo_regions!,
    regularize_field_boundary_conditions,
    FieldBoundaryConditions,
    apply_z_bcs!,
    ValueBoundaryCondition

using Oceananigans.Architectures: architecture
using Oceananigans.Utils: prettysummary, prettytime, launch!
using Oceananigans.Fields: CenterField, Field, Center, Face, interior, TracerFields
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Operators
using Oceananigans.Forcings: model_forcing

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!

struct ForwardEulerTimeStepper{T}
    tendencies :: T

    function ForwardEulerTimeStepper(grid)
        tendencies = (thickness = Field{Center, Center, Nothing}(grid),
                      enthalpy = CenterField(grid),
                      temperature = CenterField(grid),
                      salinity = CenterField(grid))
        T = typeof(tendencies)
        return new{T}(tendencies)
    end
end
    
mutable struct ThicknessCoordinateSeaIceModel{TS, G, C, H, E, S, T, F, R, Q, FT, Li} <: AbstractModel{TS}
    # Ice grid, non-dimensionalized in the vertical with z ∈ (0, 1)
    grid :: G
    clock :: C

    # Internal ice state
    ice_thickness :: H
    ice_enthalpy :: E
    ice_salinity :: S
    ice_temperature :: T

    # Equations
    forcing :: F
    penetrating_shortwave_radiation :: R
    heat_transport :: Q

    # For time-stepping
    timestepper :: TS

    # Ice properties
    ice_conductivity :: FT
    ice_density :: FT
    ice_heat_capacity :: FT
    latent_heat_of_fusion :: FT
    liquidus :: Li
end

struct HeatConduction{K}
    conductivity :: K
end

const PTTSIM = ThicknessCoordinateSeaIceModel

const ice_property_names = (:ice_conductivity,
                            :ice_density,
                            :ice_heat_capacity,
                            :latent_heat_of_fusion)

@inline gather_ice_properties(model) = NamedTuple(name => getproperty(model, name) for name in ice_property_names)

@inline ice_state(model) = (e = model.ice_enthalpy,
                            T = model.ice_temperature,
                            S = model.ice_salinity,
                            h = model.ice_thickness)

fields(model::PTTSIM) = ice_state(model)

prognostic_fields(model::PTTSIM) = (e = model.ice_enthalpy,
                                    S = model.ice_salinity,
                                    h = model.ice_thickness)

"""
    ThicknessCoordinateSeaIceModel(grid; kw...)

A model for sea ice.
"""
function ThicknessCoordinateSeaIceModel(grid;
                                        clock = Clock{eltype(grid)}(0, 0, 1),
                                        boundary_conditions = NamedTuple(),
                                        forcing = NamedTuple(),
                                        penetrating_shortwave_radiation = nothing,
                                        heat_transport = HeatConduction(2.03),
                                        timestepper = ForwardEulerTimeStepper(grid),
                                        ice_density = 911.0,
                                        ice_heat_capacity = 2000.0,
                                        ice_conductivity = 10.0, # ...
                                        latent_heat_of_fusion = 333.55e3, # J/kg
                                        liquidus = LinearLiquidus())

    # Build user-defined boundary conditions
    field_names = (:e, :T, :S, :h)
    boundary_conditions = regularize_field_boundary_conditions(boundary_conditions, grid, field_names)

    # Initialize the state
    three_dimensional_state = TracerFields((:e, :S, :T), grid, boundary_conditions)
    ice_thickness = Field{Center, Center, Nothing}(grid)
    ice_enthalpy, ice_salinity, ice_temperature = three_dimensional_state 

    prognostic_fields = (e = ice_enthalpy, S = ice_salinity, h = ice_thickness)
    forcing = model_forcing(prognostic_fields; forcing...)

    FT = eltype(grid)

    return ThicknessCoordinateSeaIceModel(grid,
                                                       clock,
                                                       ice_thickness,
                                                       ice_enthalpy,
                                                       ice_salinity,
                                                       ice_temperature,
                                                       forcing,
                                                       penetrating_shortwave_radiation,
                                                       heat_transport,
                                                       timestepper,
                                                       FT(ice_conductivity),
                                                       FT(ice_density),
                                                       FT(ice_heat_capacity),
                                                       FT(latent_heat_of_fusion),
                                                       liquidus)
end

#####
##### User interface for setting model state
#####

@inline function ice_temperature(i, j, k, grid, state, properties)
    ρi = properties.ice_density
    ci = properties.ice_heat_capacity
    e = state.e
    return @inbounds e[i, j, k] / (ρi * ci)
end

@inline function ice_enthalpy(i, j, k, grid, state, properties)
    ρi = properties.ice_density
    ci = properties.ice_heat_capacity
    T = state.T
    return @inbounds ρi * ci * T[i, j, k]
end

function compute_ice_temperature!(model)
    T = model.ice_temperature
    S = model.ice_salinity
    e = model.ice_enthalpy
    state = ice_state(model)
    properties = gather_ice_properties(model)
    grid = model.grid

    Nx, Ny, Nz = size(grid)
    @inbounds for i = 1:Nx, j = 1:Ny, k = 1:Nz
        T[i, j, k] = ice_temperature(i, j, k, grid, state, properties)
    end

    fill_halo_regions!(T, model.clock, fields(model))

    return nothing
end

function compute_ice_enthalpy!(model)
    T = model.ice_temperature
    S = model.ice_salinity
    e = model.ice_enthalpy
    properties = gather_ice_properties(model)
    state = ice_state(model)
    grid = model.grid

    Nx, Ny, Nz = size(grid)
    @inbounds for i = 1:Nx, j = 1:Ny, k = 1:Nz
        e[i, j, k] = ice_enthalpy(i, j, k, grid, state, properties)
    end

    return nothing
end

function set!(model::PTTSIM; T=nothing, e=nothing, S=nothing, h=nothing)
    isnothing(S) || set!(model.ice_salinity, S)
    isnothing(h) || set!(model.ice_thickness, h)

    # Either T or e but not both
    if !isnothing(T) && !isnothing(e)
        throw(ArgumentError("Can set either T or e, but not both!"))
    elseif !isnothing(T)
        set!(model.ice_temperature, T)
        compute_ice_enthalpy!(model)
    elseif !isnothing(e)
        set!(model.ice_enthalpy, e)
        compute_ice_temperature!(model)
    end

    update_state!(model)

    return nothing
end

#####
##### Time-stepping
#####

@inline penetrating_radiation_flux_divergence(i, j, k, grid, penetrating_shortwave_radiation::Nothing, state) = zero(grid)

@inline heat_transport_flux_divergence(i, j, k, grid, heat_transport::HeatConduction, state) =
    ∂zᶜᶜᶜ(i, j, k, grid, vertical_conductive_flux, heat_transport.conductivity, state.T)
        
@inline vertical_conductive_flux(i, j, k, grid, K::Number, T) = - K * ∂zᶜᶜᶠ(i, j, k, grid, T)

function compute_tendencies!(model)
    grid = model.grid
    arch = architecture(grid)
    K = model.ice_conductivity
    T = model.ice_temperature
    e = model.ice_enthalpy
    Ge = model.timestepper.tendencies.enthalpy
    heat_transport = model.heat_transport
    penetrating_shortwave_radiation = model.penetrating_shortwave_radiation
    forcing = model.forcing
    clock = model.clock
    properties = gather_ice_properties(model)

    model_fields = fields(model)
    state = ice_state(model)

    # Interior contribution
    Nx, Ny, Nz = size(grid)
    @inbounds for i = 1:Nx, j = 1:Ny, k = 1:Nz
        # enthalpy tendency from flux divergence
        Ge[i, j, k] = (- heat_transport_flux_divergence(i, j, k, grid, heat_transport, state)
                       + penetrating_radiation_flux_divergence(i, j, k, grid, penetrating_shortwave_radiation, state)
                       + forcing.e(i, j, k, grid, clock, model_fields))
    end

    # Calculate fluxes
    apply_z_bcs!(Ge, e, arch, model.clock, fields(model))

    return nothing
end

function step_ice_enthalpy!(model, Δt)
    grid = model.grid
    e = model.ice_enthalpy
    Ge = model.timestepper.tendencies.enthalpy

    # Forward Euler
    Nx, Ny, Nz = size(grid)
    @inbounds for i = 1:Nx, j = 1:Ny, k = 1:Nz
        e[i, j, k] = e[i, j, k] + Δt * Ge[i, j, k]
    end

    return nothing
end

function update_state!(model::PTTSIM)
    compute_ice_temperature!(model)
    return nothing
end

function time_step!(model::PTTSIM, Δt; callbacks=nothing)
    update_state!(model)
    compute_tendencies!(model)
    step_ice_enthalpy!(model, Δt)
    # step_ice_thickness!(model, Δt)
    
    tick!(model.clock, Δt)
    update_state!(model)

    return nothing
end

