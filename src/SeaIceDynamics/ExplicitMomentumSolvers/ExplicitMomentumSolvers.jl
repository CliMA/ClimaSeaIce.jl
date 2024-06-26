module ExplicitMomentumSolvers

export ExplicitMomentumSolver

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

using ClimaSeaIce.SeaIceDynamics.Rheologies: 
    compute_stresses!,
    initialize_rheology!,
    x_internal_stress_divergence,
    y_internal_stress_divergence

import ClimaSeaIce.SeaIceDynamics: 
        AbstractMomentumSolver,
        step_momentum!,
        update_stepping_coefficients!,
        get_stepping_coefficients

struct ExplicitMomentumSolver{U, R, FT, A} <: AbstractMomentumSolver
    previous_velocities :: U # ice velocities at time step n
    rheology :: R # Rheology to compute stresses
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
by substepping a `substeps` amount of time. 

# Arguments
=============
- `grid`: The grid on which the solver operates.

# Keyword Arguments
====================
- `ocean_ice_drag_coefficient`: coefficient for the ocean - ice drag, default `5.5e-3`.
- `rheology`: The rheology model used to calculate the divergence of the internal stresses ∇ ⋅ σ. Defaults to `ExplicitViscoPlasticRheology(grid)`.
- `ocean_ice_drag_coefficient`: The drag coefficient between ocean and ice. Defaults to `5.5e-3`.
- `substepping_coefficient`: Coefficient for substepping momentum (β) and internal stresses (α) (depends on the particular formulation).
              Default value is `DynamicSteppingCoefficient(grid)`.
- `substeps`: Number of substeps for the visco-plastic calculation. Default value is `400`.
              Note that we here assume that β (the modified EVP parameter that substeps velocity)
              is equal to the number of substeps.

## Returns
An instance of `ExplicitMomentumSolver` with the specified parameters.
"""
function ExplicitMomentumSolver(grid; 
                                rheology = ExplicitViscoPlasticRheology(grid),
                                ocean_ice_drag_coefficient = 5.5e-3,
                                substepping_coefficient = DynamicSteppingCoefficient(grid),
                                substeps = 400)

    uⁿ = XFaceField(grid) 
    vⁿ = YFaceField(grid) 

    previous_velocities = (u = uⁿ, v = vⁿ)

    return ExplicitMomentumSolver(previous_velocities,
                                  rheology,
                                  ocean_ice_drag_coefficient,
                                  substepping_coefficient,
                                  substeps)
end

function initialize_substepping!(model, solver::ExplicitMomentumSolver)
    uⁿ, vⁿ = solver.previous_velocities
    u,  v  = model.velocities

    launch!(architecture(model.grid), model.grid, :xy, _store_initial_velocities!, uⁿ, vⁿ, u, v)
    fill_halo_regions!((uⁿ, vⁿ), model.clock, fields(model))
    initialize_rheology!(model, solver.rheology)

    return nothing
end

# We need initial velocities for the momentum update step
@kernel function _store_initial_velocities!(uⁿ, vⁿ, u, v)
    i, j = @index(Global, NTuple)

    @inbounds uⁿ[i, j, 1] = u[i, j, 1]
    @inbounds vⁿ[i, j, 1] = v[i, j, 1]
end

include("explicit_sea_ice_dynamics.jl")
include("simple_stepping_coefficients.jl")
include("dynamic_stepping_coefficients.jl")
include("momentum_stepping_kernels.jl")

end