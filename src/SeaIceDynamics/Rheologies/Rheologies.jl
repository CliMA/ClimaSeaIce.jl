module Rheologies

export ExplicitViscoPlasticRheology

using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.Grids
using Oceananigans.Grids: architecture
using Oceananigans.Utils
using Oceananigans.Operators

using KernelAbstractions: @kernel, @index

using Adapt

using ClimaSeaIce.SeaIceDynamics

using ClimaSeaIce.SeaIceDynamics:
    AbstractRheology,
    update_stepping_coefficients!,
    get_stepping_coefficients

import Oceananigans.BoundaryConditions: fill_halo_regions!
import Oceananigans.ImmersedBoundaries: mask_immersed_field!

include("nothing_rheology.jl")
include("initialize_rheologies.jl")
include("explicit_visco_plastic_rheology.jl")
include("brittle_bingham_maxwell_rheology.jl")

end

