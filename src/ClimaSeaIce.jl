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
   
@inline ice_mass(i, j, k, grid, h, ℵ, ρ) = @inbounds h[i, j, k] * ρ[i, j, k] * ℵ[i, j, k]

# TODO: move this to Oceananigans
include("forward_euler_timestepper.jl")

include("SeaIceThermodynamics/SeaIceThermodynamics.jl")
include("Rheologies/Rheologies.jl")
include("SeaIceDynamics/SeaIceDynamics.jl")
include("sea_ice_model.jl")
include("sea_ice_advection.jl")
include("tracer_tendency_kernel_functions.jl")
include("EnthalpyMethodSeaIceModel.jl")

using .SeaIceThermodynamics
using .SeaIceDynamics
using .Rheologies

# Timestepping
include("sea_ice_fe_step.jl")
include("sea_ice_rk_substep.jl")

# Advection timescale for a `SeaIceModel`. Sea ice dynamics are two-dimensional so 
# we reuse the `cell_advection_timescale` function defined in Oceananigans by passing 
# `w = ZeroField()`.
function cell_advection_timescale(model::SeaIceModel)
    velocities = merge(model.velocities, (; w = ZeroField()))
    return cell_advection_timescale(model.grid, velocities)
end

# No diffusion timescale for sea ice for now
cell_diffusion_timescale(model::SeaIceModel) = Inf

#####
##### Default output attributes for NetCDF output
#####

default_horizontal_velocity_attributes(::RectilinearGrid) = Dict(
    "u" => Dict("long_name" => "Velocity in the +x-direction.", "units" => "m/s"),
    "v" => Dict("long_name" => "Velocity in the +y-direction.", "units" => "m/s"))

default_horizontal_velocity_attributes(::LatitudeLongitudeGrid) = Dict(
    "u" => Dict("long_name" => "Velocity in the zonal direction (+ = east).", "units" => "m/s"),
    "v" => Dict("long_name" => "Velocity in the meridional direction (+ = north).", "units" => "m/s"))

default_horizontal_velocity_attributes(::OrthogonalSphericalShellGrid) = Dict(
    "u" => Dict("long_name" => "Velocity in the i-direction (+ = increasing i).", "units" => "m/s"),
    "v" => Dict("long_name" => "Velocity in the j-direction (+ = increasing j).", "units" => "m/s"))

default_horizontal_velocity_attributes(ibg::ImmersedBoundaryGrid) = default_horizontal_velocity_attributes(ibg.underlying_grid)

default_sea_ice_attributes() = Dict(
    "h" => Dict("long_name" => "Sea ice thickness.", "units" => "m"),
    "ℵ" => Dict("long_name" => "Sea ice concentration.", "units" => "-"))

function Oceananigans.OutputWriters.default_output_attributes(model::SeaIceModel)
    velocity_attrs = default_horizontal_velocity_attributes(model.grid)
    tracer_attrs = default_sea_ice_attributes()
    return merge(velocity_attrs, tracer_attrs)
end

end # module
