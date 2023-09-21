using ClimaSeaIce:
    PhaseTransitions,
    ForwardEulerTimestepper

using ClimaSeaIce.ThermalBoundaryConditions:
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
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!, initialize!
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Utils: prettytime
import Oceananigans.Simulations: iteration

# TODO: move to Oceananigans
# import Oceananigans.Fields: field
# field(loc, a::Number, grid) = ConstantField(a)

struct SlabSeaIceModel{GR, CL, TS, IT, IC, ST, IS, U, STF, TBC, CF, P, MIT, A}
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
    external_thermal_fluxes :: STF
    thermal_boundary_conditions :: TBC
    # Internal flux
    internal_thermal_flux :: CF
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
    print(io, "└── external_thermal_fluxes: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_thermal_fluxes.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_thermal_fluxes.bottom, "     "))
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
                         clock                             = Clock{eltype(grid)}(0, 0, 1),
                         ice_thickness                     = Field{Center, Center, Nothing}(grid),
                         ice_consolidation_thickness       = 0.0, # m
                         ice_concentration                 = Field{Center, Center, Nothing}(grid),
                         ice_salinity                      = 0, # psu
                         top_surface_temperature           = nothing,
                         top_thermal_flux                  = nothing,
                         bottom_thermal_flux               = 0,
                         velocities                        = nothing,
                         advection                         = nothing,
                         top_thermal_boundary_condition    = MeltingConstrainedFluxBalance(),
                         bottom_thermal_boundary_condition = IceWaterThermalEquilibrium(),
                         internal_thermal_flux             = ConductiveFlux(eltype(grid), conductivity=2),
                         phase_transitions                 = PhaseTransitions(eltype(grid)))

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

    # Construct an internal thermal flux function that captures the liquidus and
    # bottom boundary condition.
    parameters = (flux = internal_thermal_flux,
                  liquidus = phase_transitions.liquidus,
                  bottom_thermal_boundary_condition = bottom_thermal_boundary_condition)

    internal_thermal_flux_function = FluxFunction(slab_internal_thermal_flux;
                                                  parameters,
                                                  top_temperature_dependent=true)

    # Construct default top thermal flux if one is not provided
    if isnothing(top_thermal_flux)
        if top_thermal_boundary_condition isa PrescribedTemperature  
            # Default: external top flux is in equilibrium with internal fluxes
            top_thermal_flux = internal_thermal_flux_function
        else
            # Default: no external top surface flux
            top_thermal_flux = 0
        end
    end

    # Construct default top temperature if one is not provided
    if isnothing(top_surface_temperature)
        # Check top boundary condition
        if top_thermal_boundary_condition isa PrescribedTemperature  
            top_surface_temperature = top_thermal_boundary_condition.temperature 
        else # build the default
            top_surface_temperature = Field{Center, Center, Nothing}(grid)
        end
    end

    # Convert to `field` (does nothing if it's already a Field)
    top_surface_temperature = field((Center, Center, Nothing), top_surface_temperature, grid)

    # Package the external fluxes and boundary conditions
    external_thermal_fluxes = (top = top_thermal_flux,    
                               bottom = bottom_thermal_flux) 

    thermal_boundary_conditions = (top = top_thermal_boundary_condition,
                                   bottom = bottom_thermal_boundary_condition)

    return SlabSeaIceModel(grid,
                           clock,
                           timestepper,
                           ice_thickness,
                           ice_concentration,
                           top_surface_temperature,
                           ice_salinity,
                           velocities,
                           external_thermal_fluxes,
                           thermal_boundary_conditions,
                           internal_thermal_flux_function,
                           phase_transitions,
                           ice_consolidation_thickness,
                           advection)
end

function set!(model::SSIM; h=nothing, α=nothing)
    !isnothing(h) && set!(model.ice_thickness, h)
    !isnothing(α) && set!(model.ice_conentration, α)
    return nothing
end

