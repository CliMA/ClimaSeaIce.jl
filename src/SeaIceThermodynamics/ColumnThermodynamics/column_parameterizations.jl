"""
    PrescribedBulkSalinity(profile=nothing)

Closure indicating that bulk salinity is prescribed. `profile` may be any value accepted by `set!` for an Oceananigans field.
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

function MaykutUntersteinerConductivity(FT::DataType=Oceananigans.defaults.FloatType; fresh_ice_conductivity = 2.03, salinity_coefficient = 0.13,
                                        minimum_conductivity = 0.10, temperature_floor = 1e-11)
    return MaykutUntersteinerConductivity{FT}(fresh_ice_conductivity, salinity_coefficient, minimum_conductivity, temperature_floor)
end

Base.summary(::MaykutUntersteinerConductivity) = "MaykutUntersteinerConductivity"

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

function BubblyBrineConductivity(FT::DataType=Oceananigans.defaults.FloatType; ice_density = 917, pure_ice_density = 917,
                                 minimum_conductivity = 0.10, temperature_floor = 1e-11)
    return BubblyBrineConductivity{FT}(ice_density, pure_ice_density, minimum_conductivity, temperature_floor)
end

Base.summary(::BubblyBrineConductivity) = "BubblyBrineConductivity"

@inline ice_thermal_conductivity(conductivity::Number, T, S) = conductivity

@inline function ice_thermal_conductivity(conductivity::MaykutUntersteinerConductivity, T, S)
    T⁻ = min(-conductivity.temperature_floor, T)
    K = conductivity.fresh_ice_conductivity + conductivity.salinity_coefficient * S / T⁻
    return max(K, conductivity.minimum_conductivity)
end

@inline function ice_thermal_conductivity(conductivity::BubblyBrineConductivity, T, S)
    T⁻ = min(-conductivity.temperature_floor, T)
    K = conductivity.ice_density / conductivity.pure_ice_density * (2.11 - 0.011 * T + 0.09 * S / T⁻)
    return max(K, conductivity.minimum_conductivity)
end

# Harmonic (series-resistance) mean of the two adjacent layer conductivities.
@inline function face_thermal_conductivity(conductivity, Tᵇ, Sᵇ, Δzᵇ, Tᵗ, Sᵗ, Δzᵗ)
    Kᵇ = ice_thermal_conductivity(conductivity, Tᵇ, Sᵇ)
    Kᵗ = ice_thermal_conductivity(conductivity, Tᵗ, Sᵗ)
    return (Δzᵇ + Δzᵗ) / (Δzᵇ / Kᵇ + Δzᵗ / Kᵗ)
end

"""
    ConductiveTemperatureTransport([FT=Oceananigans.defaults.FloatType; conductivity=2, diffusivity=0])

Energy transport by conduction down the temperature gradient plus direct diffusion of internal energy.
`conductivity` [W m⁻¹ K⁻¹] is a number or a conductivity closure; `diffusivity` [m² s⁻¹] multiplies the
internal-energy gradient.
"""
struct ConductiveTemperatureTransport{K, D}
    conductivity :: K
    diffusivity :: D
end

function ConductiveTemperatureTransport(FT::DataType=Oceananigans.defaults.FloatType; conductivity = 2, diffusivity = 0)
    return ConductiveTemperatureTransport(convert_parameter(FT, conductivity), convert(FT, diffusivity))
end

"""
    NoSalinityTransport()

Salinity transport closure without diffusion; prognostic bulk salinity only follows the moving grid.
"""
struct NoSalinityTransport end

"""
    BulkSalinityDiffusion([FT=Oceananigans.defaults.FloatType; diffusivity=0])

Closed-boundary diffusion of prognostic bulk salinity.
"""
struct BulkSalinityDiffusion{K}
    diffusivity :: K
end

BulkSalinityDiffusion(FT::DataType=Oceananigans.defaults.FloatType; diffusivity = 0) = BulkSalinityDiffusion(convert(FT, diffusivity))

"""
    NoShortwaveAbsorption()

Shortwave closure with no internal shortwave source.
"""
struct NoShortwaveAbsorption end

"""
    ExponentialShortwaveAbsorption([FT=Oceananigans.defaults.FloatType; surface_transmission=0, attenuation_scale=1])

Beer-law shortwave closure. `surface_transmission` is the shortwave flux transmitted through the top face
and `attenuation_scale` the vertical e-folding length.
"""
struct ExponentialShortwaveAbsorption{T, L}
    surface_transmission :: T
    attenuation_scale :: L
end

function ExponentialShortwaveAbsorption(FT::DataType=Oceananigans.defaults.FloatType; surface_transmission = 0, attenuation_scale = 1)
    return ExponentialShortwaveAbsorption(convert(FT, surface_transmission), convert(FT, attenuation_scale))
end
