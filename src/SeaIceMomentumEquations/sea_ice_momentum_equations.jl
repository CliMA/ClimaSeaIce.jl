using ClimaSeaIce.Rheologies
using Adapt

struct SeaIceMomentumEquation{S, C, R, V, A, FT}
    coriolis :: C
    rheology :: R
    auxiliary_fields :: A
    solver :: S
    free_drift_velocity :: V
    minimum_concentration :: FT
    minimum_mass :: FT
end

""" a simple explicit solver """
struct ExplicitSolver end

function SeaIceMomentumEquation(grid; 
                                coriolis = nothing,
                                rheology = ElastoViscoPlasticRheology(eltype(grid)),
                                auxiliary_fields = NamedTuple(),
                                free_drift_velocity = OceanSurfaceVelocity(grid),
                                solver = ExplicitSolver(),
                                minimum_concentration = 1e-3,
                                minimum_mass = 0.1)

    auxiliary_fields = merge(auxiliary_fields, required_auxiliary_fields(rheology, grid))

    return SeaIceMomentumEquation(coriolis, 
                                  rheology, 
                                  auxiliary_fields, 
                                  solver,
                                  free_drift_velocity,
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
