using ClimaSeaIce.Rheologies
using Adapt

struct SeaIceMomentumEquation{S, C, R, V, A, FT}
    coriolis :: C
    rheology :: R
    auxiliary_fields :: A
    solver :: S
    ocean_velocities :: V
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
                           ocean_velocities=OceanSurfaceVelocity(grid),
                           solver=ExplicitSolver(),
                           minimum_concentration=1e-3,
                           minimum_mass=1.0)

Constructs a `SeaIceMomentumEquation` object that controls the dynamical evolution of sea-ice momentum.

Arguments
=========
- `grid`: The computational grid.

Keyword Arguments
=================

- `coriolis`: Parameters for the background rotation rate of the model.
- `rheology`: The sea ice rheology model, default is `ElastoViscoPlasticRheology(eltype(grid))`.
- `auxiliary_fields`: A named tuple of auxiliary fields, default is an empty `NamedTuple()`.
- `ocean_velocities`: The ocean surface velocities used to limit the sea ice momentum when the mass or the concentration are
                      below a certain threshold. default is `OceanSurfaceVelocity(grid)`.
- `solver`: The momentum solver to be used.
- `minimum_concentration`: The minimum sea ice concentration above which the sea ice velocity is dynamically calculated, default is `1e-3`.
- `minimum_mass`: The minimum sea ice mass per area above which the sea ice velocity is dynamically calculated, default is `1.0 kg/m²`.
"""
function SeaIceMomentumEquation(grid; 
                                coriolis = nothing,
                                rheology = ElastoViscoPlasticRheology(eltype(grid)),
                                auxiliary_fields = NamedTuple(),
                                ocean_velocities = OceanSurfaceVelocity(grid),
                                solver = ExplicitSolver(),
                                minimum_concentration = 1e-3,
                                minimum_mass = 1.0)

    auxiliary_fields = merge(auxiliary_fields, required_auxiliary_fields(rheology, grid))

    return SeaIceMomentumEquation(coriolis, 
                                  rheology, 
                                  auxiliary_fields, 
                                  solver,
                                  ocean_velocities,
                                  minimum_concentration,
                                  minimum_mass)
end

fields(mom::SeaIceMomentumEquation) = mom.auxiliary_fields

struct OceanSurfaceVelocity{U, V, FT}
    u :: U
    v :: V
    C :: FT
end

""" 
    OceanSurfaceVelocity(grid; velocities=nothing, mitigation=0.01)

Build a type that controls the free drift velocity of the sea ice where concentration and mass are lower than a threshold.
The free drift velocity is calculated as `velocities` mitigated by a factor `mitigation`.
"""
function OceanSurfaceVelocity(grid; velocities=nothing, mitigation=0.01)
    if isanothing(velocities)
        u = XFaceField(grid)
        v = YFaceField(grid)
    else
        u = velocities.u
        v = velocities.v
    end

    return OceanSurfaceVelocity(u, v, mitigation)
end

@inline free_drift_u(i, j, k, grid, f::OceanSurfaceVelocity) = @inbounds f.u[i, j, k] * f.C
@inline free_drift_v(i, j, k, grid, f::OceanSurfaceVelocity) = @inbounds f.v[i, j, k] * f.C

# Just passing velocities
@inline free_drift_u(i, j, k, grid, f) = @inbounds f.u[i, j, k] 
@inline free_drift_v(i, j, k, grid, f) = @inbounds f.v[i, j, k] 
