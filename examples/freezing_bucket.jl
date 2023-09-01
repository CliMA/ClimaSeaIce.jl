using Oceananigans
using Oceananigans.Units
using ClimaSeaIce

# # Setting up a model of a bucket of water freezing from the top down
#
# We'll model a slab of ice in a bucket with one grid point, which
# requires only a zero-dimensional grid.

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

## Set the temperature at the top of the bucket to ``Tᵤ = -10ᵒC``.
model = SlabSeaIceModel(grid; top_thermal_boundary_condition=PrescribedTemperature(-10))

## We'll freeze the bucket for 10 straight days
simulation = Simulation(model, Δt=10minute, stop_time=10days)

## Initialize the ice thickness to 1 cm
set!(model, h=0.01)

# # Collecting data and running the simulation
#
# Before running the simulation, we set up a `Callback` to create
# a timeseries of the ice thickness saved at every time step.

## Container to hold the data
timeseries = []

## Callback function to collect the data from the `sim`ulation
function accumulate_timeseries(sim)
    h = sim.model.ice_thickness
    push!(timeseries, (time(sim), first(h)))
end

## Add the callback to `simulation`
simulation.callbacks[:save] = Callback(accumulate_timeseries)

# Now we're ready to hit the Big Red Button (it should run pretty quick):
run!(simulation)

# # Visualizing the result
#
# It'd be a shame to run such a "cool" simulation without looking at the
# results. We'll visualize it with Makie.

using GLMakie

# `timeseries` is a `Vector` of `Tuple`. So we have to do a bit of processing
# to build `Vector`s of time `t` and thickness `h`. It's not much work though:
t = map(ts -> ts[1], timeseries)
h = map(ts -> ts[2], timeseries)

# Just for fun, we also compute the velocity of the ice-water interface:
dhdt = @. (h[2:end] - h[1:end-1]) / simulation.Δt

# All that's left, really, is to put those `lines!` in an `Axis`:
set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(1200, 600))

axh = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Ice thickness (cm)")
axd = Axis(fig[1, 2], xlabel="Ice thickness (cm)", ylabel="Freezing rate (μm s⁻¹)")

lines!(axh, t ./ day, 1e2 .* h)
lines!(axd, 1e2 .* h[1:end-1], 1e6 .* dhdt)

display(fig)

