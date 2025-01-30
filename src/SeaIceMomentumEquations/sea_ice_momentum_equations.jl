using ClimaSeaIce.Rheologies
using Adapt

struct SeaIceMomentumEquation{S, C, R, V, A, ES, FT}
    coriolis :: C
    rheology :: R
    auxiliary_fields :: A
    solver :: S
    ocean_velocities :: V
    external_momentum_stresses :: ES
    minimum_concentration :: FT
    minimum_mass :: FT
end

""" a simple explicit solver """
struct ExplicitSolver end

"""
    SeaIceMomentumEquation(grid; 
                           coriolis=nothing,
                           rheology=ElastoViscoPlasticRheology(eltype(grid)),
                           auxiliary_fields=NamedTuple(),
                           ocean_velocities=nothing,
                           solver=ExplicitSolver(),
                           minimum_concentration=1e-3,
                           minimum_mass=1.0)

Constructs a `SeaIceMomentumEquation` object that controls the dynamical evolution of sea-ice momentum.
The sea-ice momentum obey the following evolution equation:
```math
    ∂u                   τₒ    τₐ
    -- + f x u = ∇ ⋅ σ + --  + -- 
    ∂t                   mᵢ    mᵢ
```
where the terms (left to right) represent (1) the time derivative of the ice velocity, (2) the coriolis force.
(3) the divergence of internal stresses, (4) the ice-ocean boundary stress, and (5) the ice-atmosphere boundary stress.

Arguments
=========
- `grid`: The computational grid.

Keyword Arguments
=================

- `coriolis`: Parameters for the background rotation rate of the model.
- `rheology`: The sea ice rheology model, default is `ElastoViscoPlasticRheology(eltype(grid))`.
- `auxiliary_fields`: A named tuple of auxiliary fields, default is an empty `NamedTuple()`.
- `ocean_velocities`: The ocean surface velocities used to limit the sea ice momentum when the mass or the concentration are
                      below a certain threshold. default is `nothing` (indicating that the free drift velocities are zero).
- `solver`: The momentum solver to be used.
- `minimum_concentration`: The minimum sea ice concentration above which the sea ice velocity is dynamically calculated, default is `1e-3`.
- `minimum_mass`: The minimum sea ice mass per area above which the sea ice velocity is dynamically calculated, default is `1.0 kg/m²`.
"""
function SeaIceMomentumEquation(grid; 
                                coriolis = nothing,
                                rheology = ElastoViscoPlasticRheology(eltype(grid)),
                                auxiliary_fields = NamedTuple(),
                                ocean_velocities = nothing,
                                top_momentum_stress    = (u = nothing, v = nothing),
                                bottom_momentum_stress = (u = nothing, v = nothing),
                                solver = ExplicitSolver(),
                                minimum_concentration = 1e-3,
                                minimum_mass = 1.0)

    auxiliary_fields = merge(auxiliary_fields, required_auxiliary_fields(rheology, grid))
    external_momentum_stresses = (top = top_momentum_stress,
                                  bottom = bottom_momentum_stress)

    return SeaIceMomentumEquation(coriolis, 
                                  rheology, 
                                  auxiliary_fields, 
                                  solver,
                                  ocean_velocities,
                                  external_momentum_stresses,
                                  eltype(grid)(minimum_concentration),
                                  eltype(grid)(minimum_mass))
end

fields(mom::SeaIceMomentumEquation) = mom.auxiliary_fields

struct MitigatedOceanSurfaceVelocity{U, V, FT}
    u :: U
    v :: V
    C :: FT
end

""" 
    OceanSurfaceVelocity(grid; velocities=nothing, mitigation=0.01)

Build a type that controls the free drift velocity of the sea ice where concentration and mass are lower than a threshold.
The free drift velocity is calculated as `velocities` mitigated by a factor `mitigation`.
"""
function MitigatedOceanSurfaceVelocity(; velocities=nothing, mitigation=0.01)
    u = velocities.u
    v = velocities.v
    return OceanSurfaceVelocity(u, v, mitigation)
end

@inline free_drift_u(i, j, k, grid, f::MitigatedOceanSurfaceVelocity) = @inbounds f.u[i, j, k] * f.C
@inline free_drift_v(i, j, k, grid, f::MitigatedOceanSurfaceVelocity) = @inbounds f.v[i, j, k] * f.C

# Just passing ocean velocities without mitigation
@inline free_drift_u(i, j, k, grid, f) = @inbounds f.u[i, j, k] 
@inline free_drift_v(i, j, k, grid, f) = @inbounds f.v[i, j, k] 

# Passing no velocities
@inline free_drift_u(i, j, k, grid, ::Nothing) = zero(grid)
@inline free_drift_v(i, j, k, grid, ::Nothing) = zero(grid)
