""" Ocean 🌊 Sea ice component of CliMa's Earth system model. """
module ClimaSeaIce

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: field, Field, Center, ZeroField, ConstantField
using Oceananigans.Grids: architecture
using Oceananigans.TimeSteppers: tick!, Clock, QuasiAdamsBashforth2TimeStepper, RungeKutta3TimeStepper
using Oceananigans.Utils

using KernelAbstractions: @kernel, @index

# Simulations interface
import Oceananigans: fields, prognostic_fields, prognostic_state, restore_prognostic_state!
import Oceananigans.Advection: cell_advection_timescale
import Oceananigans.Fields: set!
import Oceananigans.ImmersedBoundaries: mask_immersed_field!
import Oceananigans.Models: AbstractModel
import Oceananigans.Simulations: reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.TurbulenceClosures: cell_diffusion_timescale
import Oceananigans.Utils: prettytime

export SeaIceModel, 
       MeltingConstrainedFluxBalance,
       PrescribedTemperature,
       RadiativeEmission,
       PhaseTransitions,
       ConductiveFlux,
       FluxFunction,
       SlabSeaIceThermodynamics,
       SeaIceMomentumEquation, 
       ExplicitSolver, 
       SplitExplicitSolver, 
       SemiImplicitStress, 
       StressBalanceFreeDrift,
       ViscousRheology, 
       ElastoViscoPlasticRheology

# TODO: Move this to Oceananigans.jl
include("forward_euler_timestepper.jl")
   
@inline ice_mass(i, j, k, grid, h, ℵ, ρ) = @inbounds h[i, j, k] * ρ[i, j, k] * ℵ[i, j, k]

include("SeaIceThermodynamics/SeaIceThermodynamics.jl")
include("Rheologies/Rheologies.jl")
include("SeaIceDynamics/SeaIceDynamics.jl")
include("sea_ice_model.jl")
include("sea_ice_advection.jl")
include("tracer_tendency_kernel_functions.jl")
include("sea_ice_time_stepping.jl")
include("EnthalpyMethodSeaIceModel.jl")

using .SeaIceThermodynamics
using .SeaIceDynamics
using .Rheologies

# Advection timescale for a `SeaIceModel`. Sea ice dynamics are two-dimensional so 
# we reuse the `cell_advection_timescale` function defined in Oceananigans by passing 
# `w = ZeroField()`.
function cell_advection_timescale(model::SeaIceModel)
    velocities = merge(model.velocities, (; w = ZeroField()))
    return cell_advection_timescale(model.grid, velocities)
end

# No diffusion timescale for sea ice for now
cell_diffusion_timescale(model::SeaIceModel) = Inf

end # module
