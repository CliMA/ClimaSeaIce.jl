
"""
    PrescribedBulkSalinity(profile=nothing)

Closure indicating that bulk salinity is prescribed rather than prognosed.
`profile` may be any value accepted by `set!` for an Oceananigans field.
"""
struct PrescribedBulkSalinity{P}
    profile :: P
end

PrescribedBulkSalinity() = PrescribedBulkSalinity(nothing)

"""
    PrognosticBulkSalinity()

Closure indicating that bulk salinity is a prognostic column field.
"""
struct PrognosticBulkSalinity end

@inline convert_parameter(::Type{FT}, value::Number) where FT = convert(FT, value)
@inline convert_parameter(::Type{FT}, value) where FT = value

"""
    MaykutUntersteinerConductivity([FT=Oceananigans.defaults.FloatType; kwargs...])

CICE/Icepack BL99 `conduct = "MU71"` thermal conductivity closure.
"""
struct MaykutUntersteinerConductivity{FT}
    fresh_ice_conductivity :: FT
    salinity_coefficient :: FT
    minimum_conductivity :: FT
    temperature_floor :: FT
end

MaykutUntersteinerConductivity(FT::DataType=Oceananigans.defaults.FloatType;
                               fresh_ice_conductivity = 2.03,
                               salinity_coefficient = 0.13,
                               minimum_conductivity = 0.10,
                               temperature_floor = 1e-11) =
    MaykutUntersteinerConductivity(convert(FT, fresh_ice_conductivity),
                                   convert(FT, salinity_coefficient),
                                   convert(FT, minimum_conductivity),
                                   convert(FT, temperature_floor))

Base.summary(::MaykutUntersteinerConductivity) = "MaykutUntersteinerConductivity"

function Adapt.adapt_structure(to, conductivity::MaykutUntersteinerConductivity)
    return MaykutUntersteinerConductivity(Adapt.adapt(to, conductivity.fresh_ice_conductivity),
                                          Adapt.adapt(to, conductivity.salinity_coefficient),
                                          Adapt.adapt(to, conductivity.minimum_conductivity),
                                          Adapt.adapt(to, conductivity.temperature_floor))
end

"""
    BubblyBrineConductivity([FT=Oceananigans.defaults.FloatType; kwargs...])

CICE/Icepack BL99 `conduct = "bubbly"` thermal conductivity closure.
"""
struct BubblyBrineConductivity{FT}
    ice_density :: FT
    pure_ice_density :: FT
    minimum_conductivity :: FT
    temperature_floor :: FT
end

BubblyBrineConductivity(FT::DataType=Oceananigans.defaults.FloatType;
                        ice_density = 917,
                        pure_ice_density = 917,
                        minimum_conductivity = 0.10,
                        temperature_floor = 1e-11) =
    BubblyBrineConductivity(convert(FT, ice_density),
                            convert(FT, pure_ice_density),
                            convert(FT, minimum_conductivity),
                            convert(FT, temperature_floor))

Base.summary(::BubblyBrineConductivity) = "BubblyBrineConductivity"

function Adapt.adapt_structure(to, conductivity::BubblyBrineConductivity)
    return BubblyBrineConductivity(Adapt.adapt(to, conductivity.ice_density),
                                   Adapt.adapt(to, conductivity.pure_ice_density),
                                   Adapt.adapt(to, conductivity.minimum_conductivity),
                                   Adapt.adapt(to, conductivity.temperature_floor))
end

@inline ice_thermal_conductivity(conductivity::Number, temperature, bulk_salinity) = conductivity

@inline function ice_thermal_conductivity(conductivity::MaykutUntersteinerConductivity,
                                          temperature,
                                          bulk_salinity)
    T = min(-conductivity.temperature_floor, temperature)
    K = conductivity.fresh_ice_conductivity +
        conductivity.salinity_coefficient * bulk_salinity / T
    return max(K, conductivity.minimum_conductivity)
end

@inline function ice_thermal_conductivity(conductivity::BubblyBrineConductivity,
                                          temperature,
                                          bulk_salinity)
    T = min(-conductivity.temperature_floor, temperature)
    K = conductivity.ice_density / conductivity.pure_ice_density *
        (2.11 - 0.011 * temperature + 0.09 * bulk_salinity / T)
    return max(K, conductivity.minimum_conductivity)
end

@inline function face_thermal_conductivity(conductivity,
                                           lower_temperature,
                                           lower_salinity,
                                           lower_thickness,
                                           upper_temperature,
                                           upper_salinity,
                                           upper_thickness)
    lower_conductivity = ice_thermal_conductivity(conductivity,
                                                  lower_temperature,
                                                  lower_salinity)

    upper_conductivity = ice_thermal_conductivity(conductivity,
                                                  upper_temperature,
                                                  upper_salinity)

    return (lower_thickness + upper_thickness) /
           (lower_thickness / lower_conductivity +
            upper_thickness / upper_conductivity)
end

"""
    ConductiveTemperatureTransport([FT=Oceananigans.defaults.FloatType; conductivity=2])

Energy transport closure for conductive heat flux proportional to the
temperature gradient.
"""
struct ConductiveTemperatureTransport{K}
    conductivity :: K
end

ConductiveTemperatureTransport(FT::DataType=Oceananigans.defaults.FloatType;
                               conductivity = 2) =
    ConductiveTemperatureTransport(convert_parameter(FT, conductivity))

"""
    DiffusiveEnergyTransport([FT=Oceananigans.defaults.FloatType; diffusivity=0])

Energy transport closure for direct diffusion of internal energy.
"""
struct DiffusiveEnergyTransport{K}
    diffusivity :: K
end

DiffusiveEnergyTransport(FT::DataType=Oceananigans.defaults.FloatType;
                         diffusivity = 0) =
    DiffusiveEnergyTransport(convert(FT, diffusivity))

"""
    ConductiveAndDiffusiveEnergyTransport([FT=Oceananigans.defaults.FloatType; conductivity=2, diffusivity=0])

Energy transport closure combining conductive temperature transport and direct
internal-energy diffusion.
"""
struct ConductiveAndDiffusiveEnergyTransport{K, D}
    conductivity :: K
    diffusivity :: D
end

ConductiveAndDiffusiveEnergyTransport(FT::DataType=Oceananigans.defaults.FloatType;
                                      conductivity = 2,
                                      diffusivity = 0) =
    ConductiveAndDiffusiveEnergyTransport(convert_parameter(FT, conductivity),
                                          convert(FT, diffusivity))

"""
    NoSalinityTransport()

Salinity transport closure that leaves prognostic bulk salinity unchanged.
"""
struct NoSalinityTransport end

"""
    BulkSalinityDiffusion([FT=Oceananigans.defaults.FloatType; diffusivity=0])

Closed-boundary scalar diffusion closure for prognostic bulk salinity.
"""
struct BulkSalinityDiffusion{K}
    diffusivity :: K
end

BulkSalinityDiffusion(FT::DataType=Oceananigans.defaults.FloatType;
                      diffusivity = 0) =
    BulkSalinityDiffusion(convert(FT, diffusivity))

"""
    BrineSalinityDiffusion([FT=Oceananigans.defaults.FloatType; diffusivity=0])

Marker closure for brine-salinity diffusion. The first column implementation
uses [`BulkSalinityDiffusion`](@ref) for the scalar salinity step.
"""
struct BrineSalinityDiffusion{K}
    diffusivity :: K
end

BrineSalinityDiffusion(FT::DataType=Oceananigans.defaults.FloatType;
                       diffusivity = 0) =
    BrineSalinityDiffusion(convert(FT, diffusivity))

"""
    NoShortwaveAbsorption()

Shortwave-radiation closure with no internal shortwave source.
"""
struct NoShortwaveAbsorption end

"""
    ExponentialShortwaveAbsorption([FT=Oceananigans.defaults.FloatType; surface_transmission=0, attenuation_scale=1])

Beer-law shortwave-radiation closure. `surface_transmission` is the transmitted
shortwave flux at the top face of the column and `attenuation_scale` is the
vertical e-folding length.
"""
struct ExponentialShortwaveAbsorption{T, L}
    surface_transmission :: T
    attenuation_scale :: L
end

ExponentialShortwaveAbsorption(FT::DataType=Oceananigans.defaults.FloatType;
                               surface_transmission = 0,
                               attenuation_scale = 1) =
    ExponentialShortwaveAbsorption(convert(FT, surface_transmission),
                                   convert(FT, attenuation_scale))

