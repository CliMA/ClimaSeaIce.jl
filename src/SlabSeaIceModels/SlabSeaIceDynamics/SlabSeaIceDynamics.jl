module SlabSeaIceDynamics

using ClimaSeaIce
using ClimaSeaIce.SlabSeaIceModels

using Oceananigans
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Operators
using Oceananigans.Grids
using Oceananigans.Grids: architecture

# The only function we need to provide in SlabSeaIceDynamics.jl is the `step_momentum!` function.
import ClimaSeaIce.SlabSeaIceModels: step_momentum!

"""
    AbstractRheology

Abstract supertype for rheologies that inform the treatment of the stress divergence ∇⋅σ.
"""
abstract type AbstractRheology end


# Ice volume
@inline Vᵢ(i, j, k, grid, h, ℵ) = @inbounds h[i, j, k] * ℵ[i, j, k]

include("nothing_dynamics.jl") # nothing rheology, no sea-ice velocity
include("ExplicitRheologies/ExplicitRheologies.jl") # explicit rheologies

end