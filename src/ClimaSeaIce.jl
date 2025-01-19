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

function timestepping_coefficients(ts::RungeKutta3TimeStepper, substep) 
   if substep == 1 
      return ts.γ¹, zero(ts.γ¹)
   elseif substep == 2
      return ts.γ², ts.ζ²
   elseif substep == 3
      return ts.γ³, ts.ζ³
   end
end

function timestepping_coefficients(ts::QuasiAdamsBashforth2TimeStepper, args...) 
   χ  = ts.χ
   FT = eltype(χ)
   α  = + convert(FT, 1.5) + χ
   β  = - convert(FT, 0.5) + χ
   return α, β
end
   
@inline ice_mass(i, j, k, grid, h, ρ) = @inbounds h[i, j, k] * ρ[i, j, k]

include("SeaIceThermodynamics/SeaIceThermodynamics.jl")
include("Rheologies/Rheologies.jl")
include("SeaIceMomentumEquations/SeaIceMomentumEquations.jl")
include("sea_ice_model.jl")
include("sea_ice_advection.jl")
include("tracer_tendency_kernel_functions.jl")
include("sea_ice_time_stepping.jl")
include("sea_ice_ab2_time_stepping.jl")
include("sea_ice_rk3_time_stepping.jl")
include("EnthalpyMethodSeaIceModel.jl")

using .SeaIceThermodynamics
using .SeaIceMomentumEquations
using .Rheologies

end # module
