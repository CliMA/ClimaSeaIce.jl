
struct SlabSeaIceModel{GR, CL, TS, U, T, IT, IC, TD, D, A} <: AbstractModel{TS}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # Prognostic State
    velocities :: U
    tracers :: T
    ice_thickness :: IT
    ice_concentration :: IC
    # Thermodynamics
    sea_ice_thermodynamics :: TD
    # Dynamics
    sea_ice_dynamics :: D
    # Numerics
    advection :: A
end

