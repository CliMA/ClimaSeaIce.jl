module ThermalBoundaryConditions

"""
    PrescribedTemperature()

Thermal boundary condition indicating that the surface temperature is prescribed,
rather than diagnosed by balancing surface fluxes.

If the `surface_temperature` is below freezing, then no melting or phase change occurs
and `surface_fluxes` are ignored.

If the `surface_temperature` is at the melting temperature, then the `surface_fluxes` are used to
predict the melt rate ``dh/dt``.

(If `surface_temperature` is above the melting temperature, it is reset to the melting temperature
and a warning is thrown?)
"""
struct PrescribedTemperature{T}
    temperature :: T
end

include("bottom_thermal_boundary_conditions.jl")
include("surface_thermal_boundary_conditions.jl")
include("boundary_fluxes.jl")

end

