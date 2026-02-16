# # Perpetual night example
#
# In this example, we simulate sea ice under perpetual night conditions, where
# the ice surface emits longwave radiation but receives no incoming solar radiation.
# This example demonstrates how to:
#
#   * use `RadiativeEmission` for outgoing longwave radiation,
#   * set initial ice thickness,
#   * collect and visualize time series data.
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.
#
# ```julia
# using Pkg
# pkg"add Oceananigans, ClimaSeaIce, CairoMakie"
# ```
#
# ## The physical domain
#
# We use a single grid point to model a horizontally uniform ice slab:

using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: RadiativeEmission
using CairoMakie

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# ## Model configuration
#
# We build a model of an ice slab that has internal conductive heat fluxes
# and emits longwave radiation from its top surface. The `RadiativeEmission`
# boundary condition implements the Stefan-Boltzmann law for blackbody radiation:

ice_thermodynamics = SlabSeaIceThermodynamics(grid; top_heat_boundary_condition=MeltingConstrainedFluxBalance())

top_flux = Field{Nothing, Nothing, Nothing}(grid)
interior(top_flux) .= - 200.0
model = SeaIceModel(grid; top_heat_flux=(RadiativeEmission(), top_flux), ice_thermodynamics)

# We initialize the ice with a small thickness:

set!(model, h=0.01) # m

# ## Running a simulation
#
# We set up a simulation that runs for 40 days with a 1-hour time step:

simulation = Simulation(model, Δt=1hour, stop_time=40days)

# ## Collecting data
#
# Before running the simulation, we set up a callback to accumulate time series
# data of the ice thickness and top surface temperature:

timeseries = []

function accumulate_timeseries(sim)
    T = model.ice_thermodynamics.top_surface_temperature
    h = model.ice_thickness
    push!(timeseries, (time(sim), first(h), first(T)))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

# Now we run the simulation:

run!(simulation)

# ## Visualizing the results
#
# We extract the time series data and visualize the evolution of ice thickness
# and top surface temperature:

t = [datum[1] for datum in timeseries]
h = [datum[2] for datum in timeseries]
T = [datum[3] for datum in timeseries]

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(size=(1000, 800))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top temperature (ᵒC)")
axh = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")

lines!(axT, t ./ day, T)
lines!(axh, t ./ day, h)

current_figure() #hide

# Under perpetual night conditions, the ice will cool and grow as it loses
# heat through radiative emission at the surface.
