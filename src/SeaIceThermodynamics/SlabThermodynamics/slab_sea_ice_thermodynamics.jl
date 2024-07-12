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

struct SlabSeaIceThermodynamics{ST, HBC, CF, P, MIT} <: AbstractSeaIceThermodynamics
    top_surface_temperature :: ST
    heat_boundady_conditions :: HBC
    # Internal flux
    internal_heat_flux :: CF
    # Melting and freezing stuff
    phase_transitions :: P
    ice_consolidation_thickness :: MIT
end

const SSIM = SlabThermodynamics

Base.summary(model::SSIM) = "SlabThermodynamics"
prettytime(model::SSIM) = prettytime(model.clock.time)
iteration(model::SSIM) = model.clock.iteration

function Base.show(io::IO, model::SSIM)
    grid = model.grid
    arch = architecture(grid)
    gridname = typeof(grid).name.wrapper
    timestr = string("(time = ", prettytime(model), ", iteration = ", iteration(model), ")")

    print(io, "SlabThermodynamics{", typeof(arch), ", ", gridname, "}", timestr, '\n')
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
    SlabSeaIceThermodynamics(grid; kw...)

Pretty simple model for sea ice.
"""
function SlabThermodynamics(grid;
                            ice_consolidation_thickness    = 0.0, # m
                            top_surface_temperature        = nothing,
                            top_heat_boundary_condition    = MeltingConstrainedFluxBalance(),
                            bottom_heat_boundary_condition = IceWaterThermalEquilibrium(),
                            # Default internal flux: thermal conductivity of 2 kg m s⁻³ K⁻¹, appropriate for freshwater ice
                            internal_heat_flux             = ConductiveFlux(eltype(grid), conductivity=2),
                            phase_transitions              = PhaseTransitions(eltype(grid)))

    FT = eltype(grid)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_consolidation_thickness in a field
    ice_consolidation_thickness = field((Center, Center, Nothing), ice_consolidation_thickness, grid)

    # Construct an internal heat flux function that captures the liquidus and
    # bottom boundary condition.
    parameters = (flux = internal_heat_flux,
                  liquidus = phase_transitions.liquidus,
                  bottom_heat_boundary_condition = bottom_heat_boundary_condition)

    internal_heat_flux_function = FluxFunction(slab_internal_heat_flux;
                                               parameters,
                                               top_temperature_dependent=true)

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

    heat_boundary_conditions = (top = top_heat_boundary_condition,
                                bottom = bottom_heat_boundary_condition)

    return SlabSeaIceThermodynamics(grid,
                                    top_surface_temperature,
                                    heat_boundary_conditions,
                                    internal_heat_flux_function,
                                    phase_transitions,
                                    ice_consolidation_thickness)
end

function thermodynamically_consistent_top_heat_flux(top_heat_flux, thermodynamics :: SlabThermodynamics)   # Construct default top heat flux if one is not provided
    if isnothing(top_heat_flux)
        if thermodynamics.heat_boundary_conditions.top isa PrescribedTemperature  
            # Default: external top flux is in equilibrium with internal fluxes
            top_heat_flux = internal_heat_flux_function
        else
            # Default: no external top surface flux
            top_heat_flux = 0
        end
    end

    return top_heat_flux
end