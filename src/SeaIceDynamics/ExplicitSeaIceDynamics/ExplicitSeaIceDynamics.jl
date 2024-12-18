module ExplicitSeaIceDynamics

export ExplicitDynamics

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
    rheology_substeps,
    ∂ⱼ_σ₁ⱼ,
    ∂ⱼ_σ₂ⱼ,
    rheology_specific_forcing_x,
    rheology_specific_forcing_y
        
import ClimaSeaIce.SeaIceDynamics: step_momentum!

struct ExplicitDynamics{R, T, FT} 
    rheology :: R # Rheology to compute stresses
    auxiliary_fields :: T # auxiliary fields required for updating the velocity (like stresses, ice strength or additional velocities if on the E-grid)
    ocean_ice_drag_coefficient :: FT 
    substeps :: Int
end

"""
    ExplicitDynamics(grid; 
                           rheology = ExplicitViscoPlasticRheology(grid),
                           ocean_ice_drag_coefficient = 5.5,
                           substeps = 120)

Constructs an explicit momentum ice_dynamics for slab sea ice dynamics. The explicit ice_dynamics solves the momentum
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

Arguments
==========
- `grid`: The grid on which the ice_dynamics operates.

Keyword Arguments
==================
- `rheology`: The rheology model used to calculate the divergence of the internal stresses ∇ ⋅ σ. Defaults to `ExplicitViscoPlasticRheology(grid)`.
- `ocean_ice_drag_coefficient`: coefficient for the ocean - ice drag, it includes ocean density!!, default `5.5e-3 x 1000`.
- `substeps`: Number of substeps for the momentum ice_dynamics. Default value is `120`.
              Note that these substeps might be not be the ones that divide the time step in the stepping kernel.
              That is the role of the output of the `rheology_substeps` function (which, defaults to `substeps` in trivial rheologies).
"""
function ExplicitDynamics(grid; 
                                rheology = ExplicitViscoPlasticRheology(eltype(grid)),
                                ocean_ice_drag_coefficient = 5.5,
                                substeps = 120)

    auxiliary_fields = required_auxiliary_fields(grid, rheology)
    
    return ExplicitDynamics(rheology,
                                  auxiliary_fields,
                                  ocean_ice_drag_coefficient,
                                  substeps)
end

initialize_substepping!(model, ice_dynamics::ExplicitDynamics) = 
    initialize_rheology!(model, ice_dynamics.rheology)

include("explicit_sea_ice_dynamics.jl")
include("momentum_stepping_kernels.jl")

end