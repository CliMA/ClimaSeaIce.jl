module SeaIceDynamics

export ExplicitMomentumSolver, ExplicitViscoPlasticRheology, Slip, NoSlip,
       ℑxᴮᶜᶜᶜ,  ℑxᴮᶠᶜᶜ,  ℑxᴮᶜᶠᶜ,  ℑxᴮᶠᶠᶜ,
       ℑyᴮᶜᶜᶜ,  ℑyᴮᶠᶜᶜ,  ℑyᴮᶜᶠᶜ,  ℑyᴮᶠᶠᶜ,
       ℑxyᴮᶜᶜᶜ, ℑxyᴮᶠᶜᶜ, ℑxyᴮᶜᶠᶜ, ℑxyᴮᶠᶠᶜ,
       ∂xᶠᶜᶜ_c, ∂yᶜᶠᶜ_c, ∂xᶠᶠᶜ_c, ∂yᶠᶠᶜ_c,
       ∂xᶜᶜᶜ_U, ∂yᶜᶜᶜ_V, ∂xᶜᶠᶜ_U, ∂yᶠᶜᶜ_V 

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
# - a sea-ice density (Float)

# Ice volume
@inline ice_mass(i, j, k, grid, h, ℵ, ρ) = @inbounds h[i, j, k] * ℵ[i, j, k] * ρ

# Fallback functions common to `ExplicitMomentumSolvers` and `Rheologies`
@inline update_stepping_coefficients!(args...) = nothing
@inline get_stepping_coefficients(args...) = nothing

include("boundary_aware_operators.jl")
include("nothing_dynamics.jl") # nothing dynamics, sea-ice velocity is zero!
include("Rheologies/Rheologies.jl")
include("ExplicitMomentumSolvers/ExplicitMomentumSolvers.jl") # explicit momentum solvers

using .Rheologies
using .ExplicitMomentumSolvers

end