module ThermalBoundaryConditions

using Adapt

"""
    PrescribedTemperature()

Thermal boundary condition indicating that temperature is prescribed on the boundary.
"""
struct PrescribedTemperature{T}
    temperature :: T
end

Adapt.adapt_structure(to, pt::PrescribedTemperature) =
    PrescribedTemperature(adapt(to, pt.temperature))

include("bottom_thermal_boundary_conditions.jl")
include("top_thermal_boundary_conditions.jl")
include("boundary_fluxes.jl")

end

