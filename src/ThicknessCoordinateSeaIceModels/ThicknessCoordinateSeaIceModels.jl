module ThicknessCoordinateSeaIceModels

using Oceananigans
using Oceananigans.Units
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Models: AbstractModel
using Oceananigans.TimeSteppers: tick!

import Oceananigans.TimeSteppers: time_step!

include("slab_sea_ice_model.jl")
# include("thickness_coordinate_sea_ice_model.jl")

end # module

