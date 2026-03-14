struct ConductiveFlux{K}
    conductivity :: K
end

ConductiveFlux(FT::DataType=Oceananigans.defaults.FloatType; conductivity) = ConductiveFlux(convert(FT, conductivity))

@inline function slab_internal_heat_flux(conductive_flux::ConductiveFlux,
                                         top_surface_temperature,
                                         bottom_temperature,
                                         ice_thickness)

    k = conductive_flux.conductivity
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

struct IceSnowConductiveFlux{K, IT}
    snow_conductivity :: K
    ice_thermodynamics :: IT
end

Adapt.adapt_structure(to, f::IceSnowConductiveFlux) =
    IceSnowConductiveFlux(Adapt.adapt(to, f.snow_conductivity),
                          Adapt.adapt(to, f.ice_thermodynamics))

# Called by FluxFunction via getflux(snow.internal_heat_flux, i, j, grid, Tu, clock, fields, params)
# Computes the combined snow+ice conductive flux: F = (Tb - Tu) / (hs/ks + hi/ki)
@inline function ice_snow_conductive_flux(i, j, grid,
                                          top_surface_temperature::Number,
                                          clock, fields, parameters)
    ks = parameters.snow_conductivity
    ice_thermo = parameters.ice_thermodynamics
    ice_params = ice_thermo.internal_heat_flux.parameters

    ki = ice_params.flux.conductivity
    bottom_bc = ice_params.bottom_heat_boundary_condition
    liquidus = ice_params.liquidus

    Tu = top_surface_temperature
    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    @inbounds hi = fields.h[i, j, 1]
    @inbounds hs = fields.hs[i, j, 1]

    R = hs / ks + hi / ki
    return ifelse(R ≤ 0, zero(R), (Tb - Tu) / R)
end

# Compute interface temperature Tsi from surface temperature Tu
@inline function interface_temperature(i, j, grid, ice_thermo, ks, Tu, fields)
    ice_params = ice_thermo.internal_heat_flux.parameters
    ki = ice_params.flux.conductivity
    bottom_bc = ice_params.bottom_heat_boundary_condition
    liquidus = ice_params.liquidus

    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    @inbounds hi = fields.h[i, j, 1]
    @inbounds hs = fields.hs[i, j, 1]

    Ri = hi / ki
    Rs = hs / ks
    R  = Rs + Ri

    # Tsi = Tb + (Tu - Tb) * Ri / R = (ks/hs·Tu + ki/hi·Tb) / (ks/hs + ki/hi)
    Tsi = ifelse(R ≤ 0, Tb, Tb + (Tu - Tb) * Ri / R)

    return Tsi
end
