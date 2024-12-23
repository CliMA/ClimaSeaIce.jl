module SeaIceMomentumEquations

using ClimaSeaIce

using Oceananigans
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Operators
using Oceananigans.Grids
using Oceananigans.Grids: architecture

# Only two functions we need to extend
import ClimaSeaIce: compute_momentum_tendencies!, step_momentum!

## A Framework to solve for the ice momentum equation, in the form:
## 
##     ∂u                   τₒ    τₐ
##     -- + f x u = ∇ ⋅ σ + --  + -- + g∇η
##     ∂t                   mᵢ    mᵢ
## 
## where the terms (left to right) represent the 
## - time derivative of the ice velocity
## - coriolis force
## - divergence of internal stresses
## - ice-ocean boundary stress 
## - ice-atmosphere boundary stress 
## - ocean dynamic surface

include("nothing_dynamics.jl")
include("sea_ice_momentum_equations.jl")
include("momentum_tendencies_kernel_functions.jl")
include("compute_momentum_tendencies.jl")
include("step_momentum_equations.jl")

end