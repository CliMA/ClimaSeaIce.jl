module SeaIceMomentumEquations

# The only functions provided by the module
export compute_momentum_tendencies!, step_momentum!
export SeaIceMomentumEquation, ExplicitSolver, SplitExplicitSolver

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
                              required_auxiliary_fields, 
                              compute_stresses!,
                              initialize_rheology!,
                              compute_time_stepᶠᶜᶜ,
                              compute_time_stepᶜᶠᶜ,
                              sum_of_forcing_u,
                              sum_of_forcing_v

import Oceananigans: fields

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
step_momentum!(model, dynamics, Δt, stage) = nothing
compute_momentum_tendencies!(model, dynamics, Δt) = nothing

include("sea_ice_momentum_equations.jl")
include("sea_ice_external_stress.jl")
include("momentum_tendencies_kernel_functions.jl")
include("explicit_momentum_equations.jl")
include("split_explicit_momentum_equations.jl")

end