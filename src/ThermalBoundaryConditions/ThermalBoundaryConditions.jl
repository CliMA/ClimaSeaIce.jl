module ThermalBoundaryConditions

"""
    PrescribedTemperature()

Thermal boundary condition indicating that temperature is prescribed on the boundary.
"""
struct PrescribedTemperature{T}
    temperature :: T
end

include("bottom_thermal_boundary_conditions.jl")
include("top_thermal_boundary_conditions.jl")
include("boundary_fluxes.jl")

end

