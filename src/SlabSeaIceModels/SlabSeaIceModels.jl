module SlabSeaIceModels

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

using ClimaSeaIce.SeaIceDynamics: step_momentum!

include("slab_heat_and_tracer_fluxes.jl")
include("thickness_advection.jl")
include("slab_tendency_kernel_functions.jl")
include("slab_sea_ice_model.jl")
include("slab_time_stepping.jl")

using .SlabSeaIceDynamics

end # module

