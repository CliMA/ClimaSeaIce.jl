""" Ocean 🌊 Sea ice component of CliMa's Earth system model. """
module ClimaSeaIce

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Utils: prettysummary
using Oceananigans.TimeSteppers: Clock
using Oceananigans.Fields: field, Field, Center, ZeroField, ConstantField
using Oceananigans.TimeSteppers: tick!, QuasiAdamsBashforth2TimeStepper, RungeKutta3TimeStepper

using KernelAbstractions: @kernel, @index

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.Models: AbstractModel
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Utils: prettytime
import Oceananigans.ImmersedBoundaries: mask_immersed_field!

export SeaIceModel, 
       MeltingConstrainedFluxBalance,
       PrescribedTemperature,
       RadiativeEmission,
       PhaseTransitions,
       ConductiveFlux,
       FluxFunction,
       SlabSeaIceThermodynamics

struct ForwardEulerTimestepper end

import Oceananigans.ImmersedBoundaries: mask_immersed_field!

# TODO: move to Oceananigans
mask_immersed_field!(::ConstantField) = nothing
mask_immersed_field!(::ZeroField)     = nothing

include("SeaIceThermodynamics/SeaIceThermodynamics.jl")
include("SeaIceDynamics/SeaIceDynamics.jl")
include("sea_ice_model.jl")
include("sea_ice_advection.jl")
include("tracer_tendency_kernel_functions.jl")
include("sea_ice_time_stepping.jl")
include("sea_ice_ab2_time_stepping.jl")
include("sea_ice_rk3_time_stepping.jl")
include("EnthalpyMethodSeaIceModel.jl")

using .SeaIceThermodynamics

end # module
