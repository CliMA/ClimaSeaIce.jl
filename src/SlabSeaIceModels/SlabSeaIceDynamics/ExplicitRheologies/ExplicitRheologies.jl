module ExplicitRheologies

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Units

using KernelAbstractions: @kernel, @index
using Adapt

## A Framework to solve for the ice momentum equation explicitly, in the form:
## 
##     ∂u                   τₒ    τₐ
##     -- + f x u = ∇ ⋅ σ + --  + -- + g∇η
##     ∂t                   mᵢ    mᵢ
## 
## where the terms (left to right) represent the 
## - time derivative of the ice velocity
## - coriolis force
## - divergence of internal stresses
## - ice-ocean boundary stress (calculated in step_momentum!)
## - ice-atmosphere boundary stress (provided as an external flux)
## - ocean dynamic surface

# The only function we need to provide in SlabSeaIceDynamics.jl is the `step_momentum!` function.
import ClimaSeaIce.SlabSeaIceModels: step_momentum!

include("explicit_sea_ice_dynamics.jl")
include("modified_evp_stepping_coefficients.jl")
include("cavitating_flow_rheology.jl")
include("elasto_visco_plastic_rheology.jl")
include("momentum_stepping_kernels.jl")

end