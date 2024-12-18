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
using ClimaSeaIce.SeaIceDynamics: ice_mass

using Oceananigans.BoundaryConditions: fill_halo_regions!

include("nothing_rheology.jl")
include("explicit_visco_plastic_rheology.jl")

end

