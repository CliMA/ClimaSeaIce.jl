#####
##### Column-thermodynamics heat boundary condition marker. Behavior (Dirichlet, surface melt balance, direct
##### flux) lives here; the forcing value is read from `model.external_heat_fluxes` via `getflux`, exactly as the
##### slab does. The solve's dispatch methods live with the column solver.
#####

"""
    FluxBoundary()

Column heat boundary condition that injects the paired `external_heat_fluxes` flux directly across the face, with
no implicit surface solve and no melt cap. A resting column with `external_heat_fluxes` of zero is insulating.
"""
struct FluxBoundary end
