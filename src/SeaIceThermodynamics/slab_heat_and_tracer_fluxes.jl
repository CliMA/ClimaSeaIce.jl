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

struct IceSnowConductiveFlux{K}
    snow_conductivity :: K
    ice_conductivity :: K
end

Adapt.adapt_structure(to, f::IceSnowConductiveFlux) =
    IceSnowConductiveFlux(Adapt.adapt(to, f.snow_conductivity),
                          Adapt.adapt(to, f.ice_conductivity))

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
    ki = flux.ice_conductivity
    Tu = top_surface_temperature
    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    @inbounds hi = fields.h[i, j, 1]
    @inbounds hs = fields.hs[i, j, 1]

    R = hs / ks + hi / ki
    return ifelse(R ≤ 0, zero(R), (Tb - Tu) / R)
end

# Compute interface temperature Tsi from surface temperature Tu
# using the snow+ice resistance ratio: Tsi = Tb + (Tu - Tb) * Ri / (Rs + Ri)
@inline function interface_temperature(i, j, grid, flux::IceSnowConductiveFlux,
                                       bottom_bc, liquidus, Tu, fields)
    ki = flux.ice_conductivity
    ks = flux.snow_conductivity
    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    @inbounds hi = fields.h[i, j, 1]
    @inbounds hs = fields.hs[i, j, 1]

    Ri = hi / ki
    Rs = hs / ks
    R  = Rs + Ri

    Tsi = ifelse(R ≤ 0, Tb, Tb + (Tu - Tb) * Ri / R)

    return Tsi
end
