module SlabThermodynamics

export SlabSeaIceThermodynamics

using ClimaSeaIce.SeaIceThermodynamics: AbstractSeaIceThermodynamics

include("slab_heat_and_tracer_fluxes.jl")
include("slab_sea_ice_thermodynamics.jl")

end # module

