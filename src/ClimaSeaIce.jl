module ClimaSeaIce

import Oceananigans.TimeSteppers: time_step!

export
    MeltingConstrainedFluxBalance,
    PrescribedTemperature,
    RadiativeEmission,
    ConductiveFlux,
    FluxFunction,
    SlabSeaIceModel

# Use celsius, so
const reference_temperature = 273.15
const default_constant_sea_ice_salinity = 7 # psu

struct LinearLiquidus{FT}
    slope :: FT
end

LinearLiquidus(FT::DataType=Float64; slope=0.054) = LinearLiquidus(convert(FT, slope))

@inline melting_temperature(liquidus::LinearLiquidus, salinity) = - liquidus.slope * salinity

struct ForwardEulerTimestepper end

include("ThermalBoundaryConditions/ThermalBoundaryConditions.jl")

using .ThermalBoundaryConditions:
    IceWaterThermalEquilibrium,
    MeltingConstrainedFluxBalance,
    RadiativeEmission,
    FluxFunction,
    PrescribedTemperature

include("EnthalpyMethodSeaIceModels.jl")
include("SlabSeaIceModels.jl")

using .EnthalpyMethodSeaIceModels: EnthalpyMethodSeaIceModel
using .SlabSeaIceModels: SlabSeaIceModel, ConductiveFlux

end # module
