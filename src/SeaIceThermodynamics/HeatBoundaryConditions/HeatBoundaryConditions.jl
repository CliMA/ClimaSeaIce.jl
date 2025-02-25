module HeatBoundaryConditions

export MeltingConstrainedFluxBalance,
       PrescribedTemperature,
       RadiativeEmission,
       ConductiveFlux,
       FluxFunction

using Adapt

"""
    PrescribedTemperature()

heat boundary condition indicating that temperature is prescribed on the boundary.
"""
struct PrescribedTemperature{T}
    temperature :: T
end

Adapt.adapt_structure(to, pt::PrescribedTemperature) =
    PrescribedTemperature(adapt(to, pt.temperature))


@inline get_tracer(i, j, k, grid, T::AbstractArray) = @inbounds T[i, j, k]
@inline get_tracer(i, j, k, grid, T::Number) = T

include("bottom_heat_boundary_conditions.jl")
include("top_heat_boundary_conditions.jl")
include("boundary_fluxes.jl")

end
