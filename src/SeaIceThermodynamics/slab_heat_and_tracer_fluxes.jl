"""
    effective_conductivity_factor(thickness_categories)

Return the factor by which conduction through a slab of *mean* thickness underestimates conduction
through a sub-grid distribution of thicknesses.

Following Fichefet and Morales Maqueda (1997), the ice and snow within a cell are taken to be uniformly
distributed between zero and twice their mean, represented by `N = thickness_categories` equal-area
sub-categories of thickness ``(2i-1) h / N``. Those preserve the mean thickness, and because conduction
goes as ``1/h`` their mean flux exceeds the flux at the mean thickness by

∑ᴺᵢ₌₁ 1/(2i-1)

`N = 1` is conduction through the mean thickness. `N = 5`, the value used by LIM and SI3, gives 1.79.
Snow and ice are scaled by the same factor within a sub-category, so their series resistance scales
with it too and the correction is a single multiplicative constant on the conductivity.
"""
@inline effective_conductivity_factor(thickness_categories) = sum(1 / (2i - 1) for i in 1:thickness_categories)

struct ThicknessDependentConductivity{FT}
    minimum_shape :: FT
    maximum_shape :: FT
    transition_thickness :: FT
end

"""
    ThicknessDependentConductivity(FT = Oceananigans.defaults.FloatType;
                                   minimum_shape = 2.5, maximum_shape = 10, transition_thickness = 1)

Sub-grid thickness-distribution correction to the conductivity. The shape parameter of the assumed
thickness distribution decays from `maximum_shape` for vanishing thickness to `minimum_shape` for thick
ice over the scale `transition_thickness`, and conduction through the distribution exceeds conduction
through the mean thickness by ``s / (s - 1)``.
"""
ThicknessDependentConductivity(FT::DataType = Oceananigans.defaults.FloatType;
                               minimum_shape = 2.5, maximum_shape = 10, transition_thickness = 1) =
    ThicknessDependentConductivity(convert(FT, minimum_shape), convert(FT, maximum_shape),
                                   convert(FT, transition_thickness))

@inline itd_factor(::Nothing, h) = one(h)

@inline function itd_factor(c::ThicknessDependentConductivity, h)
    s = c.minimum_shape + (c.maximum_shape - c.minimum_shape) * exp(-h / c.transition_thickness)
    return s / (s - 1)
end

struct ConductiveFlux{K, S}
    conductivity :: K
    itd_shape :: S
end

"""
    ConductiveFlux(FT = Oceananigans.defaults.FloatType; conductivity, thickness_categories = 1, itd_shape = nothing)

Fourier conduction through the slab. `conductivity` is the material conductivity of the medium;
`thickness_categories` applies the sub-grid correction of [`effective_conductivity_factor`](@ref) to it,
so the stored `conductivity` is the *effective* one and every consumer of this struct inherits the
correction without further plumbing. `itd_shape` optionally applies the thickness-dependent correction of
[`ThicknessDependentConductivity`](@ref) at every evaluation.
"""
ConductiveFlux(FT::DataType=Oceananigans.defaults.FloatType; conductivity, thickness_categories=1, itd_shape=nothing) =
    ConductiveFlux(convert(FT, conductivity * effective_conductivity_factor(thickness_categories)), itd_shape)

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
    liquidus = parameters.liquidus
    Tu = top_surface_temperature
    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
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

IceSnowConductiveFlux(snow_conductivity, ice_conductivity) =
    IceSnowConductiveFlux(snow_conductivity, ice_conductivity, nothing)

Adapt.adapt_structure(to, f::IceSnowConductiveFlux) =
    IceSnowConductiveFlux(Adapt.adapt(to, f.snow_conductivity),
                          Adapt.adapt(to, f.ice_conductivity),
                          Adapt.adapt(to, f.itd_shape))

# Combined snow+ice conductive flux using resistors in series:
# F = (Tb - Tu) / (hs/ks + hi/ki)
# Uses the same parameter structure as slab_internal_heat_flux:
# parameters = (flux = IceSnowConductiveFlux, liquidus, bottom_heat_boundary_condition)
@inline function ice_snow_conductive_flux(i, j, grid,
                                          top_surface_temperature::Number,
                                          clock, fields, parameters)
    flux = parameters.flux
    bottom_bc = parameters.bottom_heat_boundary_condition
    liquidus = parameters.liquidus

    ks = flux.snow_conductivity
    Tu = top_surface_temperature
    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    @inbounds hi = fields.h[i, j, 1]
    @inbounds hs = fields.hs[i, j, 1]
    ki = flux.ice_conductivity * itd_factor(flux.itd_shape, hi)

    R = hs / ks + hi / ki
    return ifelse(R ≤ 0, zero(R), (Tb - Tu) / R)
end

# Compute interface temperature Tsi from surface temperature Tu
# using the snow+ice resistance ratio: Tsi = Tb + (Tu - Tb) * Ri / (Rs + Ri)
@inline function interface_temperature(i, j, grid, flux::IceSnowConductiveFlux,
                                       bottom_bc, liquidus, Tu, fields)
    ks = flux.snow_conductivity
    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    @inbounds hi = fields.h[i, j, 1]
    @inbounds hs = fields.hs[i, j, 1]
    ki = flux.ice_conductivity * itd_factor(flux.itd_shape, hi)

    Ri = hi / ki
    Rs = hs / ks
    R  = Rs + Ri

    Tsi = ifelse(R ≤ 0, Tb, Tb + (Tu - Tb) * Ri / R)

    return Tsi
end
