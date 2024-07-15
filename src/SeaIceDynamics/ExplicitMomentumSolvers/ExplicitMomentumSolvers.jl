module ExplicitMomentumSolvers

export ExplicitMomentumSolver, EGridDynamics, CGridDynamics

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

using ClimaSeaIce.SeaIceDynamics.Rheologies: 
    required_auxiliary_fields,
    ExplicitViscoPlasticRheology,
    compute_stresses!,
    initialize_rheology!,
    x_internal_stress_divergence,
    y_internal_stress_divergence,
    rheology_specific_numerical_terms_x,
    rheology_specific_numerical_terms_y,
    fill_stresses_halo_regions!

import ClimaSeaIce.SeaIceDynamics: 
        AbstractMomentumSolver,
        step_momentum!,
        update_stepping_coefficients!,
        get_stepping_coefficients

struct ExplicitMomentumSolver{G, R, T, FT, A} <: AbstractMomentumSolver{G}
    rheology :: R # Rheology to compute stresses
    auxiliary_fields :: T # auxiliary fields required for updating the velocity (like stresses or additional velocities if on the E-grid)
    ocean_ice_drag_coefficient :: FT 
    substepping_coefficient :: A
    substeps :: Int

    ExplicitMomentumSolver{G}(r::R, a::T, o::FT, c::A, s::Int) where {G, R, T, FT, A} = new{G, R, T, FT, A}(r, a, o, c, s)
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
                                dynamics_grid = CGridDynamics(),
                                rheology = ExplicitViscoPlasticRheology(eltype(grid)),
                                ocean_ice_drag_coefficient = 5.5e-3,
                                substepping_coefficient = DynamicSteppingCoefficient(grid),
                                substeps = 150)

    auxiliary_fields = required_auxiliary_fields(rheology, dynamics_grid)
    
    return ExplicitMomentumSolver{typeof(dynamics_grid)}(rheology,
                                                         auxiliary_fields,
                                                         ocean_ice_drag_coefficient,
                                                         substepping_coefficient,
                                                         substeps)
end

initialize_substepping!(model, solver::ExplicitMomentumSolver) = 
    initialize_rheology!(model, solver.rheology)

include("explicit_sea_ice_dynamics.jl")
include("simple_stepping_coefficients.jl")
include("dynamic_stepping_coefficients.jl")
include("cgrid_momentum_stepping_kernels.jl")
include("egrid_momentum_stepping_kernels.jl")

fill_velocities_halo_regions!(model, ::ExplicitMomentumSolver{<:CGridDynamics}, args...) = 
    fill_halo_regions!(model.velocities, args...)

function fill_velocities_halo_regions!(model, solver::ExplicitMomentumSolver{<:EGridDynamics}, args...)
    fill_halo_regions!(model.velocities, args...)

    û = solver.auxiliary_fields.v̂
    v̂ = solver.auxiliary_fields.û
    
    fill_halo_regions((û, v̂), args...)

    return nothing
end

end