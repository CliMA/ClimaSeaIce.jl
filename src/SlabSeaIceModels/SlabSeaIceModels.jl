module SlabSeaIceModels

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

# using RootSolvers: find_zero

include("difference_of_arrays.jl")
include("slab_heat_and_tracer_fluxes.jl")
include("slab_tendency_kernel_functions.jl")
include("slab_sea_ice_model.jl")
include("slab_time_stepping.jl")

end # module

