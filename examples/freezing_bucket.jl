# # Freezing bucket example
#
# A common laboratory experiment freezes an insulated bucket of water from the top
# down, using a metal lid to keep the top of the bucket at some constant, very cold
# temperature. In this example, we simulate such a scenario using the `SeaIceModel`.
# Here, the bucket is perfectly insulated and infinitely deep: if the `Simulation`
# is run for longer, the ice will keep freezing, and freezing, and will never run
# out of water. Also, the water in the infinite bucket is (somehow) all at the same
# temperature, in equilibrium with the ice-water interface (and therefore fixed at
# the melting temperature). This example demonstrates
#
#   * How to use `SlabSeaIceThermodynamics` with prescribed boundary conditions.
#   * How to configure internal heat conduction.
#   * How to use `FluxFunction` for frazil ice formation.
#   * How to collect and visualize time series data.
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
# ## Using `ClimaSeaIce.jl`
#
# Write

using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using CairoMakie

# to load ClimaSeaIce functions and objects into our script.
#
# ## An infinitely deep bucket with a single grid point
#
# Perhaps surprisingly, we need just one grid point to model a possibly infinitely
# thick slab of ice with `SlabSeaIceModel`. We would only need more than 1 grid point
# if our boundary conditions vary in the horizontal direction:

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# ## Model configuration
#
# We build our model of an ice slab freezing into a bucket. We start by defining
# a constant internal `ConductiveFlux` with ice conductivity:

conductivity = 2 # W m⁻¹ K⁻¹
internal_heat_flux = ConductiveFlux(; conductivity)

# Note that other units besides Celsius _can_ be used, but that requires setting
# `model.phase_transitions` with appropriate parameters. We set the ice heat capacity
# and density as well:

ice_heat_capacity = 2100 # J kg⁻¹ K⁻¹
ice_density = 900 # kg m⁻³
phase_transitions = PhaseTransitions(; ice_heat_capacity, ice_density)

# We set the top ice temperature:

top_temperature = -10 # °C
top_heat_boundary_condition = PrescribedTemperature(top_temperature)

# ## Constructing the thermodynamics
#
# We construct the ice thermodynamics using a simple slab sea ice representation:

ice_thermodynamics = SlabSeaIceThermodynamics(grid;
                                              internal_heat_flux,
                                              phase_transitions,
                                              top_heat_boundary_condition)

# ## Frazil ice formation
#
# We also prescribe a frazil ice heat flux that stops when the ice has reached a
# concentration of 1. This heat flux represents the initial ice formation from a
# liquid bucket:

@inline frazil_ice_formation(i, j, grid, Tuᵢ, clock, fields) = - (1 - fields.ℵ[i, j, 1]) # W m⁻²

bottom_heat_flux = FluxFunction(frazil_ice_formation)

# ## Building the model
#
# Then we assemble it all into a model:

model = SeaIceModel(grid; ice_thermodynamics, bottom_heat_flux)

# Note that the default bottom heat boundary condition for `SlabSeaIceThermodynamics`
# is `IceWaterThermalEquilibrium` with freshwater. That's what we want!

model.ice_thermodynamics.heat_boundary_conditions.bottom

# ## Running a simulation
#
# We're ready to freeze the bucket for 10 straight days. The ice will start forming
# suddenly due to the frazil ice heat flux and then eventually grow more slowly:

simulation = Simulation(model, Δt=10minute, stop_time=10days)

# ## Collecting data
#
# Before simulating the freezing bucket, we set up a `Callback` to create a time
# series of the ice thickness saved at every time step:

timeseries = []

function accumulate_timeseries(sim)
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    push!(timeseries, (time(sim), first(h), first(ℵ)))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

# Now we're ready to run the simulation:

run!(simulation)

# ## Visualizing the results
#
# It'd be a shame to run such a "cool" simulation without looking at the results.
# We'll visualize it with Makie:

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
axV = Axis(fig[1, 3], xlabel="Ice volume (cm)", ylabel="Freezing rate (μm s⁻¹)")

lines!(axh, t ./ day, 1e2 .* h)
lines!(axℵ, t ./ day, ℵ)
lines!(axV, 1e2 .* V[1:end-1], 1e6 .* dVdt)

save("freezing_bucket.png", fig)
nothing # hide

# ![](freezing_bucket.png)
#
# If you want more ice, you can increase `simulation.stop_time` and
# `run!(simulation)` again (or just re-run the whole script).
