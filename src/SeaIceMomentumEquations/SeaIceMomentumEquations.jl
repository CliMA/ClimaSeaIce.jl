module SeaIceMomentumEquations

# The only functions provided by the module
export compute_momentum_tendencies!, step_momentum!
export SeaIceMomentumEquation, ExplicitSolver, SplitExplicitSolver, SemiImplicitStress

using ClimaSeaIce

using Oceananigans
using KernelAbstractions: @kernel, @index
using Oceananigans.Utils: launch!
using Oceananigans.Operators
using Oceananigans.Grids
using Oceananigans.Grids: architecture

using ClimaSeaIce: ice_mass
using ClimaSeaIce.Rheologies: ∂ⱼ_σ₁ⱼ, 
                              ∂ⱼ_σ₂ⱼ, 
                              immersed_∂ⱼ_σ₁ⱼ,
                              immersed_∂ⱼ_σ₂ⱼ,
                              rheology_auxiliary_fields, 
                              compute_stresses!,
                              initialize_rheology!,
                              compute_substep_Δtᶠᶠᶜ,
                              sum_of_forcing_u,
                              sum_of_forcing_v

import Oceananigans: fields
import ClimaSeaIce.Rheologies: rheology_prognostic_tracers

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

# Fallbacks for `nothing` ice dynamics
step_momentum!(model, dynamics, Δt) = nothing
compute_momentum_tendencies!(model, dynamics, Δt) = nothing

include("sea_ice_momentum_equations.jl")
include("sea_ice_external_stress.jl")
include("momentum_tendencies_kernel_functions.jl")
include("explicit_momentum_equations.jl")
include("split_explicit_momentum_equations.jl")

end