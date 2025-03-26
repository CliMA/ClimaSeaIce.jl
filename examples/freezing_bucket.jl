# # A freezing bucket
#
# A common laboratory experiment freezes an insultated bucket of water
# from the top down, using a metal lid to keep the top of the bucket
# at some constant, very cold temperature. In this example, we simulate such
# a scenario using the `SeaIceModel`. Here, the bucket is perfectly insulated
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

# Next, we build our model of an ice slab freezing into a bucket.
# We start by defining a constant internal `ConductiveFlux` with 
# ice_conductivity 

conductivity = 2 # kg m s⁻³ K⁻¹
internal_heat_flux = ConductiveFlux(; conductivity)

# Note that other units besides Celsius _can_ be used, but that requires
# setting model.phase_transitions` with appropriate parameters.
# We set the ice heat capacity and density as well,

ice_heat_capacity = 2100 # J kg⁻¹ K⁻¹
ice_density = 900 # kg m⁻³
phase_transitions = PhaseTransitions(; ice_heat_capacity, ice_density)

# We set the top ice temperature,

top_temperature = -10 # ᵒC
top_heat_boundary_condition = PrescribedTemperature(-10)

# Construct the ice_thermodynamics of sea ice, for this we use a simple
# slab sea ice representation of ice_thermodynamics

ice_thermodynamics = SlabSeaIceThermodynamics(grid;
                                              internal_heat_flux,
                                              phase_transitions,
                                              top_heat_boundary_condition)

# We also prescribe a frazil ice heat flux that stops when the ice has reached a concentration of 1.
# This heat flux represents the initial ice formation from a liquid bucket.

@inline frazil_ice_formation(i, j, grid, Tuᵢ, clock, fields) = - (1 - fields.ℵ[i, j, 1]) # W m⁻²

bottom_heat_flux = FluxFunction(frazil_ice_formation)

# Then we assemble it all into a model.

model = SeaIceModel(grid; ice_thermodynamics, bottom_heat_flux)

# Note that the default bottom heat boundary condition for `SlabSeaIceThermodynamics` is
# `IceWaterThermalEquilibrium` with freshwater. That's what we want!

model.thermodynamics.heat_boundary_conditions.bottom

# Ok, we're ready to freeze the bucket for 10 straight days.
# The ice will start forming suddenly due to the frazil ice heat flux and then eventually
# grow more slowly.

simulation = Simulation(model, Δt=10minute, stop_time=10days)

# # Collecting data and running the simulation
#
# Before simulating the freezing bucket, we set up a `Callback` to create
# a timeseries of the ice thickness saved at every time step.

## Container to hold the data
timeseries = []

## Callback function to collect the data from the `sim`ulation
function accumulate_timeseries(sim)
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    push!(timeseries, (time(sim), first(h), first(ℵ)))
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
ℵ = [datum[3] for datum in timeseries]
V = h .* ℵ

# Just for fun, we also compute the velocity of the ice-water interface:
dVdt = @. (h[2:end] .* ℵ[2:end] - h[1:end-1] .* ℵ[1:end-1]) / simulation.Δt

# All that's left, really, is to put those `lines!` in an `Axis`:
set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(size=(1600, 700))

axh = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Ice thickness (cm)")
axℵ = Axis(fig[1, 2], xlabel="Time (days)", ylabel="Ice concentration (-)")
axV = Axis(fig[1, 3], xlabel="Ice Volume (cm)", ylabel="Freezing rate (μm s⁻¹)")

lines!(axh, t ./ day, 1e2 .* h)
lines!(axℵ, t ./ day, ℵ)
lines!(axV, 1e2 .* V[1:end-1], 1e6 .* dVdt)

save("freezing_bucket.png", fig)
nothing # hide

# ![](freezing_bucket.png)

# If you want more ice, you can increase `simulation.stop_time` and
# `run!(simulation)` again (or just re-run the whole script).

