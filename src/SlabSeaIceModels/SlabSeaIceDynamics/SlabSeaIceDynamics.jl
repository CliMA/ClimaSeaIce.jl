module SlabSeaIceDynamics

using ClimaSeaIce
using ClimaSeaIce.SlabSeaIceModels

using Oceananigans
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Operators
using Oceananigans.Grids
using Oceananigans.Grids: architecture

import ClimaSeaIce.SlabSeaIceModels: step_momentum!

"""
    AbstractRheology

Abstract supertype for rheologies that inform the treatment of the stress divergence ∇⋅σ.
"""
abstract type AbstractRheology end

include("nothing_dynamics.jl")
include("ExplicitRheologies/ExplicitRheologies.jl")

end