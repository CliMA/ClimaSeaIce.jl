module SlabSeaIceDynamics

export ExplicitMomentumSolver

using ClimaSeaIce
using ClimaSeaIce.SlabSeaIceModels

using Oceananigans
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Operators
using Oceananigans.Grids
using Oceananigans.Grids: architecture

using Adapt

# The only function we need to provide in SlabSeaIceDynamics.jl is the `step_momentum!` function.
import ClimaSeaIce.SlabSeaIceModels: step_momentum!

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