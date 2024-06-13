using ClimaSeaIce:
    PhaseTransitions,
    ForwardEulerTimestepper

using ClimaSeaIce.HeatBoundaryConditions:
    MeltingConstrainedFluxBalance,
    IceWaterThermalEquilibrium,
    PrescribedTemperature,
    FluxFunction,
    flux_summary

using Oceananigans.Utils: prettysummary
using Oceananigans.TimeSteppers: Clock
using Oceananigans.Fields: field, Field, Center, ZeroField, ConstantField

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.Models: AbstractModel
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Utils: prettytime

# TODO: move to Oceananigans
# import Oceananigans.Fields: field
# field(loc, a::Number, grid) = ConstantField(a)

struct SlabSeaIceModel{GR, CL, TS, IT, IC, ST, IS, U, STF, TBC, CF, P, MIT, A} <: AbstractModel{TS}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # State
    ice_thickness :: IT
    ice_concentration :: IC
    top_surface_temperature :: ST
    ice_salinity :: IS
    velocities :: U
    # Boundary conditions
    external_heat_fluxes :: STF
    heat_boundary_conditions :: TBC
    # Internal flux
    internal_heat_flux :: CF
    # Melting and freezing stuff
    phase_transitions :: P
    ice_consolidation_thickness :: MIT
    # Numerics
    advection :: A
end

const SSIM = SlabSeaIceModel

Base.summary(model::SSIM) = "SlabSeaIceModel"
prettytime(model::SSIM) = prettytime(model.clock.time)
iteration(model::SSIM) = model.clock.iteration

function Base.show(io::IO, model::SSIM)
    grid = model.grid
    arch = architecture(grid)
    gridname = typeof(grid).name.wrapper
    timestr = string("(time = ", prettytime(model), ", iteration = ", iteration(model), ")")

    print(io, "SlabSeaIceModel{", typeof(arch), ", ", gridname, "}", timestr, '\n')
    print(io, "├── grid: ", summary(model.grid), '\n')
    print(io, "├── top_surface_temperature: ", summary(model.top_surface_temperature), '\n')
    print(io, "├── minimium_ice_thickness: ", prettysummary(model.ice_consolidation_thickness), '\n')
    print(io, "└── external_heat_fluxes: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_heat_fluxes.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_heat_fluxes.bottom, "     "))
end
         
reset!(::SSIM) = nothing
initialize!(::SSIM) = nothing
default_included_properties(::SSIM) = tuple(:grid)

fields(model::SSIM) = (h = model.ice_thickness,
                       ℵ = model.ice_concentration,
                       Tᵤ = model.top_surface_temperature,
                       Sᵢ = model.ice_salinity)

# TODO: make this correct
prognostic_fields(model::SSIM) = fields(model)

"""
    SlabSeaIceModel(grid; kw...)

Pretty simple model for sea ice.
"""
function SlabSeaIceModel(grid;
                         clock                          = Clock{eltype(grid)}(time = 0),
                         ice_thickness                  = Field{Center, Center, Nothing}(grid),
                         ice_consolidation_thickness    = 0.0, # m
                         ice_concentration              = Field{Center, Center, Nothing}(grid),
                         ice_salinity                   = 0, # psu
                         top_surface_temperature        = nothing,
                         top_heat_flux                  = nothing,
                         bottom_heat_flux               = 0,
                         velocities                     = nothing,
                         advection                      = nothing,
                         top_heat_boundary_condition    = MeltingConstrainedFluxBalance(),
                         bottom_heat_boundary_condition = IceWaterThermalEquilibrium(),
                         # Default internal flux: thermal conductivity of 2 kg m s⁻³ K⁻¹, appropriate for freshwater ice
                         internal_heat_flux             = ConductiveFlux(eltype(grid), conductivity=2),
                         phase_transitions              = PhaseTransitions(eltype(grid)))

    if isnothing(velocities)
        velocities = (u = ZeroField(), v=ZeroField(), w=ZeroField())
    end

    # Only one time-stepper is supported currently
    timestepper = ForwardEulerTimestepper()
    FT = eltype(grid)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_salinity in a field
    ice_salinity = field((Center, Center, Nothing), ice_salinity, grid)
    ice_consolidation_thickness = field((Center, Center, Nothing), ice_consolidation_thickness, grid)

    # Construct an internal heat flux function that captures the liquidus and
    # bottom boundary condition.
    parameters = (flux = internal_heat_flux,
                  liquidus = phase_transitions.liquidus,
                  bottom_heat_boundary_condition = bottom_heat_boundary_condition)

    internal_heat_flux_function = FluxFunction(slab_internal_heat_flux;
                                               parameters,
                                               top_temperature_dependent=true)

    # Construct default top heat flux if one is not provided
    if isnothing(top_heat_flux)
        if top_heat_boundary_condition isa PrescribedTemperature  
            # Default: external top flux is in equilibrium with internal fluxes
            top_heat_flux = internal_heat_flux_function
        else
            # Default: no external top surface flux
            top_heat_flux = 0
        end
    end

    # Construct default top temperature if one is not provided
    if isnothing(top_surface_temperature)
        # Check top boundary condition
        if top_heat_boundary_condition isa PrescribedTemperature  
            top_surface_temperature = top_heat_boundary_condition.temperature 
        else # build the default
            top_surface_temperature = Field{Center, Center, Nothing}(grid)
        end
    end

    # Convert to `field` (does nothing if it's already a Field)
    top_surface_temperature = field((Center, Center, Nothing), top_surface_temperature, grid)

    # Package the external fluxes and boundary conditions
    external_heat_fluxes = (top = top_heat_flux,    
                            bottom = bottom_heat_flux) 

    heat_boundary_conditions = (top = top_heat_boundary_condition,
                                bottom = bottom_heat_boundary_condition)

    return SlabSeaIceModel(grid,
                           clock,
                           timestepper,
                           ice_thickness,
                           ice_concentration,
                           top_surface_temperature,
                           ice_salinity,
                           velocities,
                           external_heat_fluxes,
                           heat_boundary_conditions,
                           internal_heat_flux_function,
                           phase_transitions,
                           ice_consolidation_thickness,
                           advection)
end

function set!(model::SSIM; h=nothing, α=nothing)
    !isnothing(h) && set!(model.ice_thickness, h)
    !isnothing(α) && set!(model.ice_conentration, α)
    return nothing
end

