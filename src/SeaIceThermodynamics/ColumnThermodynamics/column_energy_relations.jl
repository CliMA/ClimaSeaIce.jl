"""
    AbstractColumnEnergyRelation

Supertype for the column enthalpy ↔ temperature ↔ bulk-salinity relations. Every thermodynamic map
(`internal_energy`, `temperature`, the derivatives, `liquid_fraction`, `brine_salinity`, `complete_melt_energy`)
is written once against the volumetric-coefficient accessors. A concrete relation only needs `phase_transitions`
plus its own `liquid_volumetric_heat_capacity` and `reference_volumetric_latent_heat` (which differ only in the
density that scales the liquid properties).
"""
abstract type AbstractColumnEnergyRelation end

@inline liquidus_slope(relation::AbstractColumnEnergyRelation) = relation.phase_transitions.liquidus.slope
@inline reference_temperature(relation::AbstractColumnEnergyRelation) = relation.phase_transitions.reference_temperature

@inline function solid_volumetric_heat_capacity(relation::AbstractColumnEnergyRelation)
    phase_transitions = relation.phase_transitions
    return phase_transitions.density * phase_transitions.heat_capacity
end

"""
    internal_energy(relation::AbstractColumnEnergyRelation, temperature, bulk_salinity)

Return the volumetric internal energy of a mushy sea-ice mixture with temperature `temperature`
and bulk salinity `bulk_salinity`, assuming a linear liquidus and vanishing solid-ice salinity.
"""
@inline function internal_energy(relation::AbstractColumnEnergyRelation, temperature, bulk_salinity)
    T₀ = reference_temperature(relation)
    m  = liquidus_slope(relation)
    a  = solid_volumetric_heat_capacity(relation)
    b  = liquid_volumetric_heat_capacity(relation) - a
    c  = reference_volumetric_latent_heat(relation)
    ΔT = temperature - T₀
    mS = m * bulk_salinity
    return ifelse(mS == 0, a * ΔT, a * ΔT - b * mS - c * mS / ΔT)
end

"""
    temperature(relation::AbstractColumnEnergyRelation, internal_energy, bulk_salinity)

Recover temperature from volumetric `internal_energy` and bulk salinity by taking the cold/mushy
root of the quadratic liquidus-energy relation.
"""
@inline function temperature(relation::AbstractColumnEnergyRelation, internal_energy, bulk_salinity)
    T₀ = reference_temperature(relation)
    m  = liquidus_slope(relation)
    a  = solid_volumetric_heat_capacity(relation)
    b  = liquid_volumetric_heat_capacity(relation) - a
    c  = reference_volumetric_latent_heat(relation)

    mS = m * bulk_salinity

    E_plus_bmS = internal_energy + b * mS
    discriminant = E_plus_bmS^2 + 4 * a * c * mS
    sqrt_discriminant = sqrt(discriminant)

    ΔT₁ = -2 * c * mS / (sqrt_discriminant + E_plus_bmS)
    ΔT₂ = (E_plus_bmS - sqrt_discriminant) / (2 * a)
    ΔT  = ifelse(E_plus_bmS > 0, ΔT₁, ΔT₂)

    return ifelse(mS == 0, T₀ + internal_energy / a, T₀ + ΔT)
end

@inline function temperature_denominator(relation::AbstractColumnEnergyRelation, internal_energy, bulk_salinity)
    T₀ = reference_temperature(relation)
    m  = liquidus_slope(relation)
    a  = solid_volumetric_heat_capacity(relation)
    c  = reference_volumetric_latent_heat(relation)
    mS = m * bulk_salinity
    ΔT = temperature(relation, internal_energy, bulk_salinity) - T₀
    return ifelse(mS == 0, a, a + c * mS / ΔT^2)
end

"""
    temperature_energy_derivative(relation, internal_energy, bulk_salinity)

Return `dT / dE` at fixed bulk salinity.
"""
@inline temperature_energy_derivative(relation::AbstractColumnEnergyRelation, internal_energy, bulk_salinity) =
    inv(temperature_denominator(relation, internal_energy, bulk_salinity))

"""
    temperature_salinity_derivative(relation, internal_energy, bulk_salinity)

Return `dT / dS` at fixed internal energy.
"""
@inline function temperature_salinity_derivative(relation::AbstractColumnEnergyRelation, internal_energy, bulk_salinity)
    T₀ = reference_temperature(relation)
    m  = liquidus_slope(relation)
    a  = solid_volumetric_heat_capacity(relation)
    b  = liquid_volumetric_heat_capacity(relation) - a
    c  = reference_volumetric_latent_heat(relation)
    ΔT = temperature(relation, internal_energy, bulk_salinity) - T₀
    denominator = a + c * m * bulk_salinity / ΔT^2
    return m * (b + c / ΔT) / denominator
end

"""
    liquid_fraction(relation, internal_energy, bulk_salinity)

Return the liquid fraction implied by the zero-solid-salinity closure.
"""
@inline function liquid_fraction(relation::AbstractColumnEnergyRelation, internal_energy, bulk_salinity)
    T₀ = reference_temperature(relation)
    m  = liquidus_slope(relation)
    T  = temperature(relation, internal_energy, bulk_salinity)
    return ifelse(bulk_salinity == 0, zero(T), m * bulk_salinity / (T₀ - T))
end

"""
    brine_salinity(relation, internal_energy, bulk_salinity)

Return the brine salinity implied by the linear liquidus.
"""
@inline function brine_salinity(relation::AbstractColumnEnergyRelation, internal_energy, bulk_salinity)
    T₀ = reference_temperature(relation)
    m  = liquidus_slope(relation)
    T  = temperature(relation, internal_energy, bulk_salinity)
    return (T₀ - T) / m
end

"""
    complete_melt_energy(relation, bulk_salinity)

Return the volumetric energy threshold above which sea ice with bulk salinity `bulk_salinity` is completely liquid.
"""
@inline function complete_melt_energy(relation::AbstractColumnEnergyRelation, bulk_salinity)
    m = liquidus_slope(relation)
    return reference_volumetric_latent_heat(relation) - liquid_volumetric_heat_capacity(relation) * m * bulk_salinity
end

"""
    QuadraticLiquidusEnergyRelation([FT=Oceananigans.defaults.FloatType; phase_transitions])

Thermodynamic relation for a zero-solid-salinity sea-ice mixture with a linear liquidus.
The relation maps between temperature, bulk salinity, and volumetric internal energy and
provides derivatives for semi-implicit column solves. Uses the liquid density for the liquid
volumetric properties.
"""
struct QuadraticLiquidusEnergyRelation{PT} <: AbstractColumnEnergyRelation
    phase_transitions :: PT
end

function QuadraticLiquidusEnergyRelation(FT::DataType=Oceananigans.defaults.FloatType;
                                         phase_transitions = PhaseTransitions(FT))
    return QuadraticLiquidusEnergyRelation(phase_transitions)
end

Base.summary(::QuadraticLiquidusEnergyRelation) = "QuadraticLiquidusEnergyRelation"

function Base.show(io::IO, relation::QuadraticLiquidusEnergyRelation)
    print(io, "QuadraticLiquidusEnergyRelation", '\n')
    print(io, "`-- phase_transitions: ", summary(relation.phase_transitions))
end

@inline function liquid_volumetric_heat_capacity(relation::QuadraticLiquidusEnergyRelation)
    phase_transitions = relation.phase_transitions
    return phase_transitions.liquid_density * phase_transitions.liquid_heat_capacity
end

@inline function reference_volumetric_latent_heat(relation::QuadraticLiquidusEnergyRelation)
    phase_transitions = relation.phase_transitions
    return phase_transitions.liquid_density * phase_transitions.reference_latent_heat
end

"""
    FixedDrainedIceSalinityProfile([FT=Oceananigans.defaults.FloatType; maximum_salinity=3.2,
                                    shape_parameter_a=0.407, shape_parameter_b=0.573])

CICE/Icepack BL99 fixed drained-ice salinity profile. The callable form expects normalized height above the ice base,
`zeta`, so that normalized downward depth from the surface is `s = 1 - zeta`.
"""
struct FixedDrainedIceSalinityProfile{FT}
    maximum_salinity :: FT
    shape_parameter_a :: FT
    shape_parameter_b :: FT
end

function FixedDrainedIceSalinityProfile(FT::DataType=Oceananigans.defaults.FloatType;
                                        maximum_salinity = 3.2,
                                        shape_parameter_a = 0.407,
                                        shape_parameter_b = 0.573)
    return FixedDrainedIceSalinityProfile(convert(FT, maximum_salinity),
                                          convert(FT, shape_parameter_a),
                                          convert(FT, shape_parameter_b))
end

Base.summary(::FixedDrainedIceSalinityProfile) = "FixedDrainedIceSalinityProfile"

function Base.show(io::IO, profile::FixedDrainedIceSalinityProfile{FT}) where FT
    print(io, "FixedDrainedIceSalinityProfile{", FT, "}", '\n')
    print(io, "|-- maximum_salinity: ", profile.maximum_salinity, '\n')
    print(io, "|-- shape_parameter_a: ", profile.shape_parameter_a, '\n')
    print(io, "`-- shape_parameter_b: ", profile.shape_parameter_b)
end

function Adapt.adapt_structure(to, profile::FixedDrainedIceSalinityProfile)
    return FixedDrainedIceSalinityProfile(Adapt.adapt(to, profile.maximum_salinity),
                                          Adapt.adapt(to, profile.shape_parameter_a),
                                          Adapt.adapt(to, profile.shape_parameter_b))
end

@inline function salinity_at_normalized_depth(profile::FixedDrainedIceSalinityProfile,
                                              normalized_depth_from_surface)
    s = normalized_depth_from_surface
    exponent = profile.shape_parameter_a / (s + profile.shape_parameter_b)
    return profile.maximum_salinity / 2 *
           (one(s) - cos(oftype(s, pi) * s^exponent))
end

@inline salinity_at_normalized_height(profile::FixedDrainedIceSalinityProfile, normalized_height_above_base) = 
    salinity_at_normalized_depth(profile, one(normalized_height_above_base) - normalized_height_above_base)

@inline (profile::FixedDrainedIceSalinityProfile)(normalized_height_above_base) =
    salinity_at_normalized_height(profile, normalized_height_above_base)

"""
    FixedSalinityBrinePocketEnergyRelation([FT=Oceananigans.defaults.FloatType; phase_transitions])

Fixed-salinity drained-ice thermodynamic relation for reproducing CICE/Icepack BL99 column physics. 
The shifted internal energy is `E = rho_i c_i T - rho_i (c_w - c_i) mu S - rho_i L_0 mu S / T`.
"""
struct FixedSalinityBrinePocketEnergyRelation{PT} <: AbstractColumnEnergyRelation
    phase_transitions :: PT
end

function FixedSalinityBrinePocketEnergyRelation(FT::DataType=Oceananigans.defaults.FloatType;
                                                phase_transitions = PhaseTransitions(FT;
                                                    heat_capacity = 2106,
                                                    liquid_heat_capacity = 4218))
    return FixedSalinityBrinePocketEnergyRelation(phase_transitions)
end

Base.summary(::FixedSalinityBrinePocketEnergyRelation) = "FixedSalinityBrinePocketEnergyRelation"

function Base.show(io::IO, relation::FixedSalinityBrinePocketEnergyRelation)
    print(io, "FixedSalinityBrinePocketEnergyRelation", '\n')
    print(io, "`-- phase_transitions: ", summary(relation.phase_transitions))
end

# BL99 scales the liquid volumetric properties by the ice density rather than the liquid density.
@inline function liquid_volumetric_heat_capacity(relation::FixedSalinityBrinePocketEnergyRelation)
    phase_transitions = relation.phase_transitions
    return phase_transitions.density * phase_transitions.liquid_heat_capacity
end

@inline function reference_volumetric_latent_heat(relation::FixedSalinityBrinePocketEnergyRelation)
    phase_transitions = relation.phase_transitions
    return phase_transitions.density * phase_transitions.reference_latent_heat
end
