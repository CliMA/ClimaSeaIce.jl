module SeaIceDynamics

export ExplicitMomentumSolver

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
# - an 
"""
    AbstractMomentumSolver

Abstract supertype for solvers of the momentum equation. Could be explicit or impicit.
For example an explicit solver resort to substepping the momentum equation within a larger
tracer advection step.
"""
abstract type AbstractMomentumSolver end

"""
    AbstractRheology

Abstract supertype for rheologies that calculate the stress divergence ∇⋅σ. 
"""
abstract type AbstractRheology end

# Ice volume
@inline Vᵢ(i, j, k, grid, h, ℵ) = @inbounds h[i, j, k] * ℵ[i, j, k]

# Fallback functions common to `ExplicitMomentumSolvers` and `Rheologies`
update_stepping_coefficients!(args...) = nothing
get_stepping_coefficients(args...) = nothing

include("nothing_dynamics.jl") # nothing dynamics, sea-ice velocity is zero!
include("Rheologies/Rheologies.jl")
include("ExplicitMomentumSolvers/ExplicitMomentumSolvers.jl") # explicit momentum solvers

end