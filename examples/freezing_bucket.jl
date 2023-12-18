# # A freezing bucket
#
# A common laboratory experiment freezes an insultated bucket of water
# from the top down, using a metal lid to keep the top of the bucket
# at some constant, very cold temperature. In this example, we simulate such
# a scenario using `SlabSeaIceModel`. Here, the bucket is perfectly insulated
# and infinitely deep, like many buckets are: if the `Simulation` is run for longer,
# the ice will keep freezing, and freezing, and will never run out of water.
# Also, the water in the infinite bucket is (somehow) all at the same temperature,
# in equilibrium with the ice-water interface (and therefore fixed at the melting
# temperature). Yep, this kind of thing happens _all the time_.
#
# We start by `using Oceananigans` to bring in functions for building grids
# and `Simulation`s and the like.

using Oceananigans
using Oceananigans.Units

# Next we `using ClimaSeaIce` to get some ice-specific names.

using ClimaSeaIce

# # An infinitely deep bucket with a single grid point
#
# Perhaps surprisingly, we need just one grid point
# to model an possibly infinitely thick slab of ice with `SlabSeaIceModel`.
# We would only need more than 1 grid point if our boundary conditions
# vary in the horizontal direction.

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# Then our model of an ice slab freezing into a bucket. We set the top ice temperature
# to `-10 ᵒC`. Note that other units besides Celsius _can_ be used, but that requires
# setting model.phase_transitions` with appropriate parameters.

model = SlabSeaIceModel(grid; top_heat_boundary_condition=PrescribedTemperature(-10))

# The default bottom heat boundary condition for `SlabSeaIceModel` is
# `IceWaterThermalEquilibrium` with freshwater. That's what we want!

model.heat_boundary_conditions.bottom

# Ok, we're ready to freeze the bucket for 10 straight days with an initial ice
# thickness of 1 cm,

simulation = Simulation(model, Δt=10minute, stop_time=10days)

set!(model, h=0.01)

# # Collecting data and running the simulation
#
# Before simulating the freezing bucket, we set up a `Callback` to create
# a timeseries of the ice thickness saved at every time step.

## Container to hold the data
timeseries = []

## Callback function to collect the data from the `sim`ulation
function accumulate_timeseries(sim)
    h = sim.model.thickness
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

using CairoMakie

# `timeseries` is a `Vector` of `Tuple`. So we have to do a bit of processing
# to build `Vector`s of time `t` and thickness `h`. It's not much work though:
t = [datum[1] for datum in timeseries]
h = [datum[2] for datum in timeseries]

# Just for fun, we also compute the velocity of the ice-water interface:
dhdt = @. (h[2:end] - h[1:end-1]) / simulation.Δt

# All that's left, really, is to put those `lines!` in an `Axis`:
set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(1200, 600))

axh = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Ice thickness (cm)")
axd = Axis(fig[1, 2], xlabel="Ice thickness (cm)", ylabel="Freezing rate (μm s⁻¹)")

lines!(axh, t ./ day, 1e2 .* h)
lines!(axd, 1e2 .* h[1:end-1], 1e6 .* dhdt)

current_figure() # hide
fig

# If you want more ice, you can increase `simulation.stop_time` and
# `run!(simulation)` again (or just re-run the whole script).

