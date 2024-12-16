struct ConductiveFlux{K}
    conductivity :: K
end

ConductiveFlux(FT::DataType=Float64; conductivity) = ConductiveFlux(convert(FT, conductivity))

@inline function slab_internal_heat_flux(conductive_flux::ConductiveFlux,
                                         top_surface_temperature,
                                         bottom_temperature,
                                         ice_thickness)

    k = conductive_flux.conductivity
    Tu = top_surface_temperature
    Tb = bottom_temperature
    h = ice_thickness

    return ifelse(h <= 0, zero(h), - k * (Tu - Tb) / h)
end

@inline function slab_internal_heat_flux(i, j, grid,
                                         top_surface_temperature::Number,
                                         clock, fields, parameters)
    flux = parameters.flux
    bottom_bc = parameters.bottom_heat_boundary_condition
    liquidus = parameters.liquidus
    Tu = top_surface_temperature
    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    h = @inbounds fields.h[i, j, 1]
    return slab_internal_heat_flux(flux, Tu, Tb, h)
end

