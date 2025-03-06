""" Ocean 🌊 Sea ice component of CliMa's Earth system model. """
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

# TODO: Move this to Oceananigans.jl
include("forward_euler_timestepper.jl")
   
@inline ice_mass(i, j, k, grid, h, ℵ, ρ) = @inbounds h[i, j, k] * ρ[i, j, k] * ℵ[i, j, k]

include("SeaIceThermodynamics/SeaIceThermodynamics.jl")
include("Rheologies/Rheologies.jl")
include("SeaIceMomentumEquations/SeaIceMomentumEquations.jl")
include("sea_ice_model.jl")
include("sea_ice_advection.jl")
include("tracer_tendency_kernel_functions.jl")
include("sea_ice_time_stepping.jl")
include("EnthalpyMethodSeaIceModel.jl")

using .SeaIceThermodynamics
using .SeaIceMomentumEquations
using .Rheologies

end # module
