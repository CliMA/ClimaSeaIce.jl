struct UniformThicknessDistribution{FT}
    conductivity_factor :: FT
end

"""
    UniformThicknessDistribution(FT = Oceananigans.defaults.FloatType; categories = 5)

Sub-grid thickness distribution of Fichefet and Morales Maqueda (1997): the ice and snow within a cell are
uniformly distributed between zero and twice their mean, represented by `N = categories` equal-area
sub-categories of thickness ``(2i-1) h / N``. Those preserve the mean thickness, and because conduction
goes as ``1/h`` their mean flux exceeds the flux at the mean thickness by

∑ᴺᵢ₌₁ 1/(2i-1)

`N = 1` is conduction through the mean thickness. `N = 5`, the value used by LIM and SI3, gives 1.79.
"""
function UniformThicknessDistribution(FT::DataType = Oceananigans.defaults.FloatType; categories = 5)
    conductivity_factor = sum(1 / (2i - 1) for i in 1:categories)
    return UniformThicknessDistribution(convert(FT, conductivity_factor))
end

struct GammaThicknessDistribution{FT}
    minimum_shape :: FT
    maximum_shape :: FT
    transition_thickness :: FT
end

"""
    GammaThicknessDistribution(FT = Oceananigans.defaults.FloatType;
                               minimum_shape = 2.5, maximum_shape = 10, transition_thickness = 1)

Sub-grid thickness distribution taken as a gamma distribution whose shape parameter ``s`` decays from
`maximum_shape` for vanishing mean thickness to `minimum_shape` for thick ice over the scale
`transition_thickness`. Since conduction goes as ``1/h``, conduction through a gamma distribution exceeds
conduction through its mean thickness by ``s / (s - 1)``.
"""
function GammaThicknessDistribution(FT::DataType = Oceananigans.defaults.FloatType;
                                    minimum_shape = 2.5, maximum_shape = 10, transition_thickness = 1)
    return GammaThicknessDistribution(convert(FT, minimum_shape), convert(FT, maximum_shape), convert(FT, transition_thickness))
end

@inline itd_factor(::Nothing, h) = one(h)
@inline itd_factor(d::UniformThicknessDistribution, h) = d.conductivity_factor

@inline function itd_factor(c::GammaThicknessDistribution, h)
    s = c.minimum_shape + (c.maximum_shape - c.minimum_shape) * exp(-h / c.transition_thickness)
    return s / (s - 1)
end

struct ConductiveFlux{K, S}
    conductivity :: K
    itd_shape :: S
end

"""
    ConductiveFlux(FT = Oceananigans.defaults.FloatType; conductivity, itd_shape = nothing)

Fourier conduction through the slab with the material `conductivity`, multiplied by the `itd_factor` of the
sub-grid thickness distribution `itd_shape`: `nothing`, a `UniformThicknessDistribution`, or a
`GammaThicknessDistribution`.
"""
function ConductiveFlux(FT::DataType=Oceananigans.defaults.FloatType; conductivity, itd_shape=nothing)
    return ConductiveFlux(convert(FT, conductivity), itd_shape)
end

@inline function slab_internal_heat_flux(conductive_flux::ConductiveFlux,
                                         top_surface_temperature,
                                         bottom_temperature,
                                         ice_thickness)

    k = conductive_flux.conductivity * itd_factor(conductive_flux.itd_shape, ice_thickness)
    Tu = top_surface_temperature
    Tb = bottom_temperature
    h = ice_thickness

    return ifelse(h ≤ 0, zero(h), - k * (Tu - Tb) / h)
end

@inline function slab_internal_heat_flux(i, j, grid,
                                         top_surface_temperature::Number,
                                         clock, fields, parameters)
    flux = parameters.flux
    bottom_bc = parameters.bottom_heat_boundary_condition
    phase_transitions = parameters.phase_transitions
    Tu = top_surface_temperature
    Tb = bottom_temperature(i, j, grid, bottom_bc, phase_transitions, fields)
    hi = @inbounds fields.h[i, j, 1]
    return slab_internal_heat_flux(flux, Tu, Tb, hi)
end

#####
##### IceSnowConductiveFlux — combined resistors-in-series for the snow layer
#####

struct IceSnowConductiveFlux{K, S}
    snow_conductivity :: K
    ice_conductivity :: K
    itd_shape :: S
end

IceSnowConductiveFlux(snow_conductivity, ice_conductivity) = IceSnowConductiveFlux(snow_conductivity, ice_conductivity, nothing)

Adapt.adapt_structure(to, f::IceSnowConductiveFlux) = IceSnowConductiveFlux(Adapt.adapt(to, f.snow_conductivity),
                                                                            Adapt.adapt(to, f.ice_conductivity),
                                                                            Adapt.adapt(to, f.itd_shape))

# Combined snow+ice conductive flux using resistors in series, scaled by the sub-grid thickness factor f
# that applies to snow and ice alike: F = f (Tb - Tu) / (hs/ks + hi/ki)
# Uses the same parameter structure as slab_internal_heat_flux:
# parameters = (flux = IceSnowConductiveFlux, phase_transitions, bottom_heat_boundary_condition)
@inline function ice_snow_conductive_flux(i, j, grid,
                                          top_surface_temperature::Number,
                                          clock, fields, parameters)
    flux = parameters.flux
    bottom_bc = parameters.bottom_heat_boundary_condition
    phase_transitions = parameters.phase_transitions

    ks = flux.snow_conductivity
    ki = flux.ice_conductivity
    Tu = top_surface_temperature
    Tb = bottom_temperature(i, j, grid, bottom_bc, phase_transitions, fields)
    @inbounds hi = fields.h[i, j, 1]
    @inbounds hs = fields.hs[i, j, 1]
    f = itd_factor(flux.itd_shape, hi)

    R = hs / ks + hi / ki
    return ifelse(R ≤ 0, zero(R), f * (Tb - Tu) / R)
end

# Compute interface temperature Tsi from surface temperature Tu
# using the snow+ice resistance ratio: Tsi = Tb + (Tu - Tb) * Ri / (Rs + Ri)
@inline function interface_temperature(i, j, grid, flux::IceSnowConductiveFlux,
                                       bottom_bc, phase_transitions, Tu, fields)
    ki = flux.ice_conductivity
    ks = flux.snow_conductivity
    Tb = bottom_temperature(i, j, grid, bottom_bc, phase_transitions, fields)
    @inbounds hi = fields.h[i, j, 1]
    @inbounds hs = fields.hs[i, j, 1]

    Ri = hi / ki
    Rs = hs / ks
    R  = Rs + Ri

    Tsi = ifelse(R ≤ 0, Tb, Tb + (Tu - Tb) * Ri / R)

    return Tsi
end
