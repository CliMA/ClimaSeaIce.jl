""" Ocean 🌊 Sea ice component of CliMa's Earth system model. """
module ClimaSeaIce

using Oceananigans
using Oceananigans.Utils: prettysummary
using Oceananigans.TimeSteppers: Clock
using Oceananigans.Fields: field, Field, Center, ZeroField, ConstantField

using KernelAbstractions: @kernel, @index

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.Models: AbstractModel
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Utils: prettytime

export SeaIceModel, 
       MeltingConstrainedFluxBalance,
       PrescribedTemperature,
       RadiativeEmission,
       PhaseTransitions,
       ConductiveFlux,
       FluxFunction,
       SlabSeaIceThermodynamics

struct ForwardEulerTimestepper end

include("SeaIceThermodynamics/SeaIceThermodynamics.jl")
include("SeaIceDynamics/SeaIceDynamics.jl")
include("sea_ice_model.jl")
include("sea_ice_advection.jl")
include("tracer_tendency_kernel_functions.jl")
include("sea_ice_time_stepping.jl")
include("EnthalpyMethodSeaIceModel.jl")

using .SeaIceThermodynamics
using .SeaIceDynamics

end # module
