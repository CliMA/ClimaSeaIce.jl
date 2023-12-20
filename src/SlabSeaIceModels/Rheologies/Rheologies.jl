module Rheologies

using ClimaSeaIce
using ClimaSeaIce.SlabSeaIceModels

using Oceananigans
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Operators
using Oceananigans.Grids
using Oceananigans.Grids: architecture

import ClimaSeaIce.SlabSeaIceModels: step_momentum!
import ClimaSeaIce.SlabSeaIceModels: SlabSeaIceModelTendencyFields

## A Framework to solve for the ice momentum equation, in the form:
## 
##     ∂u                          τₒ    τₐ
##     -- + f x (u - uₒ) = ∇ ⋅ σ + --  + --
##     ∂t                          mᵢ    mᵢ
## 
## where the terms (left to right) represent the 
## - time derivative of the ice velocity
## - coriolis force
## - divergence of internal stresses
## - ice-ocean boundary stress
## - ice-atmosphere boundary stress

"""
    AbstractRheology

Abstract supertype for rheologies that inform the treatment of the stress divergence ∇⋅σ.
"""
abstract type AbstractRheology end


include("difference_of_arrays.jl")
include("momentum_stepping_kernels.jl")
include("nothing_rheology.jl")
include("explicit_rheology.jl")
include("free_drift_rheology.jl")
include("cavitating_flow_rheology.jl")
# include("elasto_visco_plastic_rheology.jl")

end