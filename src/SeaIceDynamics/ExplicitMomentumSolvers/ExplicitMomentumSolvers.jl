module ExplicitMomentumSolvers

export ExplicitMomentumSolver

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Units
using Oceananigans.BoundaryConditions
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Oceananigans.Operators

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

using ClimaSeaIce.SeaIceDynamics
using ClimaSeaIce.SeaIceDynamics: ice_mass

using ClimaSeaIce.SeaIceDynamics.Rheologies: 
    ExplicitViscoPlasticRheology,
    compute_stresses!,
    initialize_rheology!,
    required_auxiliary_fields,
    ∂ⱼ_σ₁ⱼ,
    ∂ⱼ_σ₂ⱼ,
    rheology_specific_forcing_x,
    rheology_specific_forcing_y
    
import ClimaSeaIce.SeaIceDynamics: 
        step_momentum!,
        rhoelogy_substeps

struct ExplicitMomentumSolver{R, B, T, FT, A} 
    rheology :: R # Rheology to compute stresses
    auxiliary_fields :: T # auxiliary fields required for updating the velocity (like stresses, ice strength or additional velocities if on the E-grid)
    ocean_ice_drag_coefficient :: FT 
    substepping_coefficient :: A
    substeps :: Int
end

"""
    ExplicitMomentumSolver(grid; 
                           rheology = ExplicitViscoPlasticRheology(grid),
                           ocean_ice_drag_coefficient = 5.5e-3,
                           substepping_coefficient = DynamicSteppingCoefficient(grid),
                           substeps = 1000)

Constructs an explicit momentum solver for slab sea ice dynamics. The explicit solver solves the momentum
equation for sea ice:

```math
   ∂u                   τₒ    τₐ
   -- + f x u = ∇ ⋅ σ + --  + -- + g∇η
   ∂t                   mᵢ    mᵢ
```
by substepping a `substeps` amount of time. In practice the substepping solves:

```math
  uᵖ⁺¹ - uᵖ = β⁻¹ * (Δt / mᵢ * (∇ ⋅ σᵖ⁺¹ + fk̂ × uᵖ + τₐ + τₒ) + uⁿ - uᵖ)
```
where `Δt` is the large advective time step and ``β`` is a `stepping_coefficient`
designed to obtain convergence.

# Arguments
=============
- `grid`: The grid on which the solver operates.

# Keyword Arguments
====================
- `ocean_ice_drag_coefficient`: coefficient for the ocean - ice drag, it includes ocean density!!, default `5.5e-3 x 1000`.
- `rheology`: The rheology model used to calculate the divergence of the internal stresses ∇ ⋅ σ. Defaults to `ExplicitViscoPlasticRheology(grid)`.
- `ocean_ice_drag_coefficient`: The drag coefficient between ocean and ice. Defaults to `5.5e-3`.
- `substepping_coefficient`: Coefficient for substepping momentum (β) and internal stresses (α) (depends on the particular formulation).
              Default value is `DynamicSteppingCoefficient(grid)`.
- `substeps`: Number of substeps for the momentum solver. Default value is `120`.
              Note that we here assume that β (the modified EVP parameter that substeps velocity)
              is equal to the number of substeps.

## Returns
An instance of `ExplicitMomentumSolver` with the specified parameters.
"""
function ExplicitMomentumSolver(grid; 
                                rheology = ExplicitViscoPlasticRheology(eltype(grid)),
                                ocean_ice_drag_coefficient = 5.5,
                                substepping_coefficient = DynamicSteppingCoefficient(grid),
                                substeps = 120)

    auxiliary_fields = required_auxiliary_fields(grid, rheology)
    
    return ExplicitMomentumSolver(rheology,
                                  auxiliary_fields,
                                  ocean_ice_drag_coefficient,
                                  substepping_coefficient,
                                  substeps)
end

initialize_substepping!(model, solver::ExplicitMomentumSolver) = 
    initialize_rheology!(model, solver.rheology)

include("explicit_sea_ice_dynamics.jl")
include("simple_stepping_coefficients.jl")
include("momentum_stepping_kernels.jl")

end