module SeaIceDynamics

# The only functions provided by the module
export compute_momentum_tendencies!, time_step_momentum!
export SeaIceMomentumEquation, ExplicitSolver, SplitExplicitSolver, SemiImplicitStress, StressBalanceFreeDrift

using Adapt: Adapt
using KernelAbstractions: @kernel, @index
using Oceananigans: Oceananigans, prognostic_state, prognostic_fields,
                    restore_prognostic_state!, fields
using Oceananigans.Architectures: architecture
using Oceananigans.DistributedComputations: Distributed
using Oceananigans.Fields: Field
using Oceananigans.Grids: Center, Face
using Oceananigans.Operators: ℑxyᶜᶠᵃ, ℑxyᶠᶜᵃ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ
using Oceananigans.TimeSteppers: SplitRungeKuttaTimeStepper
using Oceananigans.Utils: KernelParameters, launch!

using ..ClimaSeaIce: ice_mass
using ..Rheologies: ∂ⱼ_σ₁ⱼ, ∂ⱼ_σ₂ⱼ,
                    immersed_∂ⱼ_σ₁ⱼ, immersed_∂ⱼ_σ₂ⱼ,
                    Auxiliaries,compute_stresses!,
                    initialize_rheology!, finalize_rheology!,
                    compute_substep_Δtᶠᶜᶜ, compute_substep_Δtᶜᶠᶜ,
                    sum_of_forcing_u, sum_of_forcing_v

## A Framework to solve for the ice momentum equation, in the form:
##
##     ∂u/∂t + f x u = ∇ ⋅ σ / mᵢ  + τₒ / mᵢ + τₐ / mᵢ + g ∇η
##
## where the terms (left to right) represent the
## - time derivative of the ice velocity
## - coriolis force
## - divergence of internal stresses
## - ice-ocean boundary stress
## - ice-atmosphere boundary stress
## - ocean dynamic surface

# Fallbacks for `nothing` ice dynamics
time_step_momentum!(model, dynamics, Δt) = nothing
compute_momentum_tendencies!(model, dynamics, Δt) = nothing

include("sea_ice_external_stress.jl")
include("stress_balance_free_drift.jl")
include("sea_ice_momentum_equations.jl")
include("momentum_tendencies_kernel_functions.jl")
include("explicit_momentum_equations.jl")
include("split_explicit_momentum_equations.jl")

end
