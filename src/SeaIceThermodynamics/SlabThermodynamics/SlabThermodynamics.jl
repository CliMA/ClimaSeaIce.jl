module SlabThermodynamics

export SlabSeaIceThermodynamics

using ClimaSeaIce:
    ForwardEulerTimestepper

using ClimaSeaIce.SeaIceThermodynamics:
    PhaseTransitions
    AbstractSeaIceThermodynamics,
    latent_heat,

using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions:
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

include("slab_heat_and_tracer_fluxes.jl")
include("slab_sea_ice_thermodynamics.jl")

end # module

