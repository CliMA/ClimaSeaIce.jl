""" Ocean 🌊 Sea ice component of CliMa's Earth system model. """
module ClimaSeaIce

import Oceananigans.TimeSteppers: time_step!

export
    MeltingConstrainedFluxBalance,
    PrescribedTemperature,
    RadiativeEmission,
    PhaseTransitions,
    ConductiveFlux,
    FluxFunction,
    SlabSeaIceModel

#####
##### A bit of thermodynamics to start the day
#####

struct LinearLiquidus{FT}
    freshwater_melting_temperature :: FT
    slope :: FT
end

"""
    LinearLiquidus(FT=Float64,
                   slope = 0.054, # psu / ᵒC
                   freshwater_melting_temperature = 0) # ᵒC

Return a linear model for the dependence of the melting temperature of
saltwater on salinity,

```math
Tₘ(S) = T₀ - m S ,
```

where ``Tₘ(S)`` is the melting temperature as a function of salinity ``S``,
``T₀`` is the melting temperature of freshwater, and ``m`` is the ratio
between the melting temperature and salinity (in other words the linear model
should be thought of as defining ``m`` and could be written ``m ≡ (T₀ - Tₘ) / S``.
The signs are arranged so that ``m > 0`` for saltwater).

The defaults assume that salinity is given in practical salinity units `psu` and
temperature is in degrees Celsius.

Note: the function `melting_temperature(liquidus, salinity)` returns the
melting temperature given `salinity`.
"""
function LinearLiquidus(FT::DataType=Float64;
                        slope = 0.054, # psu / ᵒC
                        freshwater_melting_temperature = 0) # ᵒC

    return LinearLiquidus(convert(FT, freshwater_melting_temperature),
                          convert(FT, slope))
end

@inline function melting_temperature(liquidus::LinearLiquidus, salinity)
    return liquidus.freshwater_melting_temperature - liquidus.slope * salinity
end

struct PhaseTransitions{FT, L}
    ice_density :: FT
    ice_heat_capacity :: FT
    liquid_density :: FT
    liquid_heat_capacity :: FT
    reference_latent_heat :: FT
    reference_temperature :: FT
    liquidus :: L
end

"""
    PhaseTransitions(FT=Float64,
                     ice_density           = 917,   # kg m⁻³
                     ice_heat_capacity     = 2000,  # J / (kg ᵒC)
                     liquid_density        = 999.8, # kg m⁻³
                     liquid_heat_capacity  = 4186,  # J / (kg ᵒC)
                     reference_latent_heat = 334e3  # J kg⁻³
                     liquidus = LinearLiquidus(FT)) # default assumes psu, ᵒC

Return a representation of transitions between the solid and liquid phases
of salty water: in other words, the freezing and melting of sea ice.

The latent heat of fusion ``ℒ(T)`` (more simply just "latent heat") is
a function of temperature ``T`` via

```math
ρᵢ ℒ(T) = ρᵢ ℒ₀ + (ρ_ℓ c_ℓ - ρᵢ cᵢ) (T - T₀)    
```

where ``ρᵢ`` is the `ice_density`, ``ρ_ℓ`` is the liquid density,
``cᵢ`` is the heat capacity of ice, and ``c_ℓ`` is the heat capacity of
liquid, and ``T₀`` is a reference temperature, all of which are assumed constant.

The default `liquidus` assumes that salinity has practical salinity units (psu)
and that temperature is degrees Celsius.
"""
@inline function PhaseTransitions(FT=Float64;
                                  ice_density           = 917,    # kg m⁻³
                                  ice_heat_capacity     = 2000,   # J / (kg ᵒC)
                                  liquid_density        = 999.8,  # kg m⁻³
                                  liquid_heat_capacity  = 4186,   # J / (kg ᵒC)
                                  reference_latent_heat = 334e3,  # J kg⁻³
                                  reference_temperature = 0,      # ᵒC
                                  liquidus = LinearLiquidus(FT))

    return PhaseTransitions(convert(FT, ice_density),
                            convert(FT, ice_heat_capacity),
                            convert(FT, liquid_density),
                            convert(FT, liquid_heat_capacity),
                            convert(FT, reference_latent_heat),
                            convert(FT, reference_temperature),
                            liquidus)
end

@inline function latent_heat(thermo::PhaseTransitions, T)
    T₀ = thermo.reference_temperature    
    ℒ₀ = thermo.reference_latent_heat
    ρᵢ = thermo.ice_density
    ρℓ = thermo.liquid_density
    cᵢ = thermo.ice_heat_capacity
    cℓ = thermo.liquid_heat_capacity

    return ρℓ * ℒ₀ + (ρℓ * cℓ - ρᵢ * cᵢ) * (T - T₀)
end

include("HeatBoundaryConditions/HeatBoundaryConditions.jl")

using .HeatBoundaryConditions:
    IceWaterThermalEquilibrium,
    MeltingConstrainedFluxBalance,
    RadiativeEmission,
    FluxFunction,
    PrescribedTemperature

include("EnthalpyMethodSeaIceModels.jl")
include("SlabSeaIceModels/SlabSeaIceModels.jl")

using .EnthalpyMethodSeaIceModels: EnthalpyMethodSeaIceModel
using .SlabSeaIceModels: SlabSeaIceModel, ConductiveFlux

end # module
