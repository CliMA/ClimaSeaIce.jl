module IceSurfaces

using Oceananigans.Fields: ZeroField, ConstantField

#####
##### Ice-atmosphere boundary condition numerical methods
#####

# Solver modes for surface temperature
struct TimeLinearizedSurfaceTemperatureSolver end

struct PrescribedTemperatureSurface{T}
    temperature :: T
end

PrescribedTemperatureSurface(T::Number) = PrescribedTemperatureSurface(ConstantField(T))

struct MeltingRadiatingSurface{S, FT, Q}
    emissivity :: FT
    stefan_boltzmann_constant :: FT
    reference_temperature :: FT
    extrinsic_heat_flux :: Q
    temperature_solver :: S
end


"""
    MeltingRadiatingSurface(FT=Float64, kw...)

Return an abstraction representing a Stefan boundary condition that estimates
the surface temperature using a balance of fluxes across the surface
(including outgoing radiation) and under the constraint that the
surface temperature is at or below the melting temperature.
If the surface temperature is equal to the melting temperature,
this implies that a flux imbalance exists, which is used to estimate the motion
of the interface.
"""
function MeltingRadiatingSurface(FT=Float64;
                                 emissivity = one(FT),
                                 stefan_boltzmann_constant = convert(FT, 5.67e-8),
                                 reference_temperature = convert(FT, 273.15),
                                 extrinsic_heat_flux = ZeroField(),
                                 temperature_solver = TimeLinearizedSurfaceTemperatureSolver())

    return MeltingRadiatingSurface(emissivity,
                                   stefan_boltzmann_constant,
                                   reference_temperature,
                                   extrinsic_heat_flux,
                                   temperature_solver)
end

struct IceWaterSurface{S}
    water_salinity :: S
end

end # module
