""" Ocean 🌊 Sea ice component of CliMa's Earth system model. """
module ClimaSeaIce

import Oceananigans.TimeSteppers: time_step!

export
    MeltingConstrainedFluxBalance,
    PrescribedTemperature,
    RadiativeEmission,
    PhaseTransitions,
    ConductiveFlux,
    FluxFunction,
    SlabSeaIceModel

struct ForwardEulerTimestepper end

include("SeaIceThermodynamics/SeaIceThermodynamics.jl")
include("sea_ice_model.jl")
include("tracer_tendency_kernel_functions.jl")
include("time_stepping.jl")

using .SeaIceThermodynamics

end # module
