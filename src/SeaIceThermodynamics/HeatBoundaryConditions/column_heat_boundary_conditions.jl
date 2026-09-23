"""
    FluxBoundary()

Column heat boundary condition that injects the paired `external_heat_fluxes` flux directly across the face, with
no implicit surface solve and no melt cap. A resting column with `external_heat_fluxes` of zero is insulating.
"""
struct FluxBoundary end
