# # Simple freezing bucket example
#
# This is ClimaSeaIce.jl's simplest example: a single grid point model of sea ice
# freezing in a bucket. This example demonstrates
#
#   * How to load `ClimaSeaIce.jl`.
#   * How to instantiate a `SlabSeaIceModel` with prescribed boundary conditions.
#   * How to configure internal heat conduction.
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.
#
# ```julia
# using Pkg
# pkg"add Oceananigans, ClimaSeaIce"
# ```
#
# ## Using `ClimaSeaIce.jl`
#
# Write

using Oceananigans
using Oceananigans.Units
using ClimaSeaIce

# to load ClimaSeaIce functions and objects into our script.
#
# ## A single grid point model
#
# For a simple freezing bucket, we only need a single grid point since we're
# modeling a horizontally uniform ice slab. We use a `Flat` topology in all
# directions to create a zero-dimensional grid:

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# ## Model configuration
#
# We configure the model with a prescribed top temperature and internal heat
# conduction. The top temperature is kept constant at -10°C, simulating a cold
# metal lid on top of the bucket:

top_temperature = -10 # °C
top_heat_boundary_condition = PrescribedTemperature(top_temperature)

# We specify internal heat conduction through the ice with a thermal conductivity
# of 2 W m⁻¹ K⁻¹ (typical for sea ice):

conductivity = 2 # W m⁻¹ K⁻¹
internal_heat_flux = ConductiveFlux(; conductivity)

# ## Building the model
#
# We assemble these components into a `SlabSeaIceModel`:

model = SlabSeaIceModel(grid; top_heat_boundary_condition, internal_heat_flux)

# The model is now ready to use! By default, the bottom boundary condition
# is `IceWaterThermalEquilibrium`, which maintains the ice-water interface
# at the melting temperature.
