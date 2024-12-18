""" Ocean 🌊 Sea ice component of CliMA's Earth system model. """
module ClimaSeaIce

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Utils: prettysummary
using Oceananigans.TimeSteppers: Clock
using Oceananigans.Fields: field, Field, Center, ZeroField, ConstantField
using Oceananigans.TimeSteppers: tick!, QuasiAdamsBashforth2TimeStepper, RungeKutta3TimeStepper
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: architecture

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
include("Rheologies.jl")
include("Advection.jl")
include("SeaIceModels/SeaIceModels.jl")
include("EnthalpyMethodSeaIceModel.jl")

using .SeaIceThermodynamics
using .SeaIceModels

end # module

