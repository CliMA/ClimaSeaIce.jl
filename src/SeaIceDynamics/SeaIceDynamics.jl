module SeaIceDynamics

export SplitExplicitDynamics, ExplicitViscoPlasticRheology

using ClimaSeaIce

using Oceananigans
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Operators
using Oceananigans.Grids
using Oceananigans.Grids: architecture

using Adapt

# The only function we need to provide from the SeaIceDynamics.jl module
# is the `step_momentum!` function. This module assumes that any model that uses
# SeaIceDynamics has 
# - 2D velocities : `u` and `v`
# - 2D ocean velocities : `u` and `v` (at the surface)
# - 2D thickness field : `h`
# - 2D concentration field `ℵ`
# - a sea-ice density (Float)

# Ice volume
@inline ice_mass(i, j, k, grid, h, ℵ, ρ) = @inbounds h[i, j, k] * ℵ[i, j, k] * ρ

include("nothing_dynamics.jl") # nothing dynamics, sea-ice velocity is zero!
include("Rheologies/Rheologies.jl")
include("SplitExplicitSeaIceDynamics/SplitExplicitSeaIceDynamics.jl") # explicit momentum solvers

using .Rheologies
using .ExplicitDynamicss

end