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
using ClimaSeaIce.SeaIceThermodynamics: prescribed_salinity_enthalpy_thermodynamics,
                                        SeaIceColumnDiscretization,
                                        ConductiveTemperatureTransport,
                                        IceWaterThermalEquilibrium
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: RadiativeEmission
import ClimaSeaIce.SeaIceThermodynamics: initialize_column_interfaces!
using CairoMakie

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# ## Model configuration
#
# We build a model of an ice slab that has internal conductive heat fluxes
# and emits longwave radiation from its top surface. The `RadiativeEmission`
# boundary condition implements the Stefan-Boltzmann law for blackbody radiation:

ice_thermodynamics = SlabThermodynamics(grid; top_heat_boundary_condition=MeltingConstrainedFluxBalance())

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

# ## Comparison with the resolved column thermodynamics
#
# The column thermodynamics reads the *same* `model.external_heat_fluxes` as the slab, so we can drive a resolved
# column with identical forcing — the outgoing `RadiativeEmission` plus the prescribed incoming flux — and compare
# it against the slab. (Both use the upward-positive flux convention, so the same flux tuple drives them
# identically.) The column resolves the internal temperature profile, while the slab assumes a steady linear one.

column_grid = RectilinearGrid(size=16, z=SeaIceColumnDiscretization((0, 2)), topology=(Flat, Flat, Bounded))

column_thermodynamics = prescribed_salinity_enthalpy_thermodynamics(column_grid;
    salinity_profile = 0.0,
    energy_transport = ConductiveTemperatureTransport(conductivity = 2),
    heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(),
                                bottom = IceWaterThermalEquilibrium(salinity = 0)))

column_model = SeaIceModel(column_grid;
                           top_heat_flux = (RadiativeEmission(), -200.0),
                           ice_thermodynamics = column_thermodynamics)

set!(column_model, h=0.01)
initialize_column_interfaces!(column_grid, column_model.ice_thickness)
set!(column_thermodynamics; bulk_salinity=0.0, temperature=-10.0)

column_simulation = Simulation(column_model, Δt=1hour, stop_time=40days)

column_timeseries = []

function accumulate_column(sim)
    Tc = sim.model.ice_thermodynamics.auxiliary.surface_temperature
    hc = sim.model.ice_thickness
    push!(column_timeseries, (time(sim), first(hc), first(Tc)))
end

column_simulation.callbacks[:save] = Callback(accumulate_column)
run!(column_simulation)

tc = [datum[1] for datum in column_timeseries]
hc = [datum[2] for datum in column_timeseries]
Tc = [datum[3] for datum in column_timeseries]

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(size=(1000, 800))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top temperature (ᵒC)")
axh = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")

lines!(axT, t  ./ day, T,  label="slab")
lines!(axT, tc ./ day, Tc, linestyle=:dash, label="resolved column")
lines!(axh, t  ./ day, h,  label="slab")
lines!(axh, tc ./ day, hc, linestyle=:dash, label="resolved column")
axislegend(axT, position=:rt)
axislegend(axh, position=:lt)

save("perpetual_night.png", fig)
nothing # hide

# ![](perpetual_night.png)
#
# Under perpetual night conditions, the ice cools and grows as it loses heat through radiative emission at the
# surface. Both thermodynamics respond to the same forcing and follow the same qualitative evolution, but the
# resolved column grows somewhat more slowly and runs a colder surface: it stores internal energy and resolves
# the temperature gradient that the zero-dimensional slab assumes is linear over the full thickness.
