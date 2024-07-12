module SlabThermodynamics

export SlabSeaIceThermodynamics

abstract type AbstractSeaIceThermodynamics end

include("slab_heat_and_tracer_fluxes.jl")
include("slab_sea_ice_thermodynamics.jl")

end # module

