module SlabSeaIceModels

using ClimaSeaIce:
    PhaseTransitions,
    latent_heat,
    ForwardEulerTimestepper,
    melting_temperature

using ClimaSeaIce.ThermalBoundaryConditions:
    MeltingConstrainedFluxBalance,
    IceWaterThermalEquilibrium,
    PrescribedTemperature,
    FluxFunction,
    SurfaceTemperatureDependent,
    bottom_temperature,
    surface_temperature,
    surface_flux_imbalance,
    bottom_flux_imbalance

# using RootSolvers: find_zero

using Oceananigans.Fields: Field, Center, ZeroField, ConstantField
using Oceananigans.TimeSteppers: Clock, tick!

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: field, set!
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!, initialize!
import Oceananigans.OutputWriters: default_included_properties

# TODO: move to Oceananigans
field(loc, a::Number, grid) = ConstantField(a)

struct SlabSeaIceModel{GR, CL, TS, IT, ST, IS, STF, TBC, CF, P}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # State
    ice_thickness :: IT
    surface_temperature :: ST
    ice_salinity :: IS
    # Boundary conditions
    external_thermal_fluxes :: STF
    thermal_boundary_conditions :: TBC
    # Internal flux
    internal_thermal_flux :: CF
    # Melting and freezing stuff
    phase_transitions :: P
end

const SSIM = SlabSeaIceModel

initialize!(::SSIM) = nothing
default_included_properties(::SSIM) = tuple(:grid)

fields(model::SSIM) = (h = model.ice_thickness,
                       T = model.surface_temperature,
                       S = model.ice_salinity)

# TODO: make this correct
prognostic_fields(model::SSIM) = fields(model)

struct ConductiveFlux{K}
    conductivity :: K
end

ConductiveFlux(FT::DataType=Float64; conductivity) = ConductiveFlux(convert(FT, conductivity))

@inline function slab_internal_thermal_flux(conductive_flux::ConductiveFlux,
                                            surface_temperature,
                                            bottom_temperature,
                                            ice_thickness)

    k = conductive_flux.conductivity
    Ts = surface_temperature
    Tb = bottom_temperature
    h = ice_thickness

    return - k * (Ts - Tb) / h
end

@inline function slab_internal_thermal_flux(i, j, grid,
                                            surface_temperature::Number,
                                            clock, fields, parameters)
    flux = parameters.flux
    bottom_bc = parameters.bottom_thermal_boundary_condition
    liquidus = parameters.liquidus
    Ts = surface_temperature
    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    h = @inbounds fields.h[i, j, 1]
    return slab_internal_thermal_flux(flux, Ts, Tb, h)
end

"""
    SlabSeaIceModel(grid; kw...)

Pretty simple model for sea ice.
"""
function SlabSeaIceModel(grid;
                         clock                              = Clock{eltype(grid)}(0, 0, 1),
                         ice_thickness                      = Field{Center, Center, Nothing}(grid),
                         ice_salinity                       = 0, # psu
                         surface_temperature                = nothing,
                         surface_thermal_flux               = 0,
                         bottom_thermal_flux                = 0,
                         surface_thermal_boundary_condition = MeltingConstrainedFluxBalance(),
                         bottom_thermal_boundary_condition  = IceWaterThermalEquilibrium(),
                         internal_thermal_flux              = ConductiveFlux(eltype(grid), conductivity=2),
                         phase_transitions                  = PhaseTransitions(eltype(grid)))

    # Only one time-stepper is supported currently
    timestepper = ForwardEulerTimestepper()
    FT = eltype(grid)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?

    # Wrap ice_salinity in a field
    ice_salinity = field((Center, Center, Nothing), ice_salinity, grid)

    # Construct default surface temperature if one is not provided
    if isnothing(surface_temperature)
        # Check surface boundary condition
        if surface_thermal_boundary_condition isa PrescribedTemperature  
            surface_temperature = surface_thermal_boundary_condition.temperature 
        else # build the default
            surface_temperature = Field{Center, Center, Nothing}(grid)
        end
    end

    # Convert to `field` (does nothing if it's already a Field)
    surface_temperature = field((Center, Center, Nothing), surface_temperature, grid)

    # Construct an internal thermal flux function that captures the liquidus and
    # bottom boundary condition.
    parameters = (flux = internal_thermal_flux,
                  liquidus = phase_transitions.liquidus,
                  bottom_thermal_boundary_condition = bottom_thermal_boundary_condition)

    internal_thermal_flux_function = FluxFunction(slab_internal_thermal_flux;
                                                  parameters,
                                                  surface_temperature_dependent=true)

    # Package the external fluxes and boundary conditions
    external_thermal_fluxes = (surface = surface_thermal_flux,    
                               bottom = bottom_thermal_flux) 

    thermal_boundary_conditions = (surface = surface_thermal_boundary_condition,
                                   bottom = bottom_thermal_boundary_condition)

    return SlabSeaIceModel(grid,
                           clock,
                           timestepper,
                           ice_thickness,
                           surface_temperature,
                           ice_salinity,
                           external_thermal_fluxes,
                           thermal_boundary_conditions,
                           internal_thermal_flux_function,
                           phase_transitions)
end

function set!(model::SSIM; h=nothing)
    !isnothing(h) && set!(model.ice_thickness, h)
    return nothing
end

function time_step!(model::SSIM, Δt; callbacks=nothing)
    grid = model.grid
    Nx, Ny, Nz = size(grid)

    phase_transitions = model.phase_transitions
    liquidus = phase_transitions.liquidus

    h = model.ice_thickness
    Ts = model.surface_temperature
    model_fields = fields(model)
    ice_salinity = model.ice_salinity
    clock = model.clock
    bottom_thermal_bc = model.thermal_boundary_conditions.bottom
    surface_thermal_bc = model.thermal_boundary_conditions.surface

    Qi = model.internal_thermal_flux
    Qb = model.external_thermal_fluxes.bottom
    Qs = model.external_thermal_fluxes.surface

    for i = 1:Nx, j=1:Ny

        @inbounds begin
            Tsⁿ = Ts[i, j, 1]
            S = ice_salinity[i, j, 1]
        end

        # 1. Update thickness with ForwardEuler step
        # 1a. Calculate residual fluxes across the top surface and bottom
        δQs = surface_flux_imbalance(i, j, grid, surface_thermal_bc, Tsⁿ, Qi, Qs, clock, model_fields)
        δQb =  bottom_flux_imbalance(i, j, grid, bottom_thermal_bc,  Tsⁿ, Qi, Qb, clock, model_fields)

        # 1b. Impose an implicit flux balance if Ts ≤ Tm.
        Tₘ = melting_temperature(liquidus, S) 
        δQs = ifelse(Tsⁿ <= Tₘ, zero(grid), δQs)

        # 1c. Calculate the latent heat of fusion at the surface and bottom
        Tbⁿ = bottom_temperature(i, j, grid, bottom_thermal_bc, liquidus)
        ℒ_Ts = latent_heat(phase_transitions, Tsⁿ)
        ℒ_Tb = latent_heat(phase_transitions, Tbⁿ)

        # 1d. Increment the thickness
        @inbounds begin
            h⁺ = h[i, j, 1] + Δt * (δQb / ℒ_Tb + δQs / ℒ_Ts)
            h⁺ = max(zero(grid), h⁺)
            h[i, j, 1] = h⁺
        end

        if !isa(surface_thermal_bc, PrescribedTemperature)
            # 2. Update surface temperature
            Ts⁺ = surface_temperature(i, j, grid, surface_thermal_bc, Tsⁿ, Qi, Qs, clock, model_fields)

            # 3. Set surface temperature (could skip this for `PrescribedTemperature`)
            @inbounds Ts[i, j, 1] = Ts⁺ 
        end
    end

    tick!(model.clock, Δt)

    return nothing
end

function update_state!(model::SSIM)

    surface_thermal_bc = model.thermal_boundary_conditions.surface

    if !isa(surface_thermal_bc, PrescribedTemperature)
        grid = model.grid
        Nx, Ny, Nz = size(grid)

        Ts = model.surface_temperature
        model_fields = fields(model)
        clock = model.clock
        surface_thermal_bc = model.thermal_boundary_conditions.surface
        Qi = model.internal_thermal_flux
        Qs = model.external_thermal_fluxes.surface

        # Update surface temperature
        for i = 1:Nx, j=1:Ny
            Ts⁻ = @inbounds Ts[i, j, 1]
            Tsⁿ = surface_temperature(i, j, grid, surface_thermal_bc, Ts⁻, Qi, Qs, clock, model_fields)
            @inbounds Ts[i, j, 1] = Tsⁿ
        end

    end

    return nothing
end

end # module
