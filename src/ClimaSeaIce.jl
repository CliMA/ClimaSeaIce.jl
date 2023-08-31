module ClimaSeaIce

import Oceananigans.TimeSteppers: time_step!

export
    MeltingConstrainedFluxBalance,
    PrescribedTemperature,
    SlabSeaIceModel

# Use celsius, so
const reference_temperature = 273.15
const default_constant_sea_ice_salinity = 7 # psu

struct LinearLiquidus{FT}
    slope :: FT
end

LinearLiquidus(FT::DataType=Float64; slope=0.054) = LinearLiquidus(convert(FT, slope))

@inline melting_temperature(liquidus::LinearLiquidus, salinity) = - liquidus.slope * salinity

struct ConstantBulkSalinity{FT}
    constant :: FT
end

ConstantBulkSalinity(FT::DataType=Float64) = 
    ConstantBulkSalinity{FT}(convert(FT, default_constant_sea_ice_salinity))

struct ForwardEulerTimestepper end

include("ThermalBoundaryConditions/ThermalBoundaryConditions.jl")

using .ThermalBoundaryConditions:
    IceWaterThermalEquilibrium,
    MeltingConstrainedFluxBalance,
    PrescribedTemperature

include("EulerianThermodynamicSeaIceModels.jl")
include("ThicknessCoordinateSeaIceModels/ThicknessCoordinateSeaIceModels.jl")

using .EulerianThermodynamicSeaIceModels: EulerianThermodynamicSeaIceModel
using .ThicknessCoordinateSeaIceModels: SlabSeaIceModel

end # module
