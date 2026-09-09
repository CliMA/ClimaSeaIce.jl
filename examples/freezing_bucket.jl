# # Freezing bucket example
#
# A common laboratory experiment freezes an insulated bucket of water from the top
# down, using a metal lid to keep the top of the bucket at some constant, very cold
# temperature. In this example, we simulate such a scenario using the `SeaIceModel`.
# Here, the bucket is perfectly insulated and infinitely deep: if the `Simulation`
# is run for longer, the ice will keep freezing, and freezing, and will never run
# out of water. Also, the water in the infinite bucket is (somehow) all at the same
# temperature, in equilibrium with the ice-water interface (and therefore fixed at
# the melting temperature). This example demonstrates how to:
#
#   * use `SlabThermodynamics` with prescribed boundary conditions,
#   * configure internal heat conduction,
#   * use `FluxFunction` for frazil ice formation,
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
# ## Using `ClimaSeaIce.jl`
#
# Write

using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: prescribed_salinity_enthalpy_thermodynamics,
                                        SeaIceColumnDiscretization,
                                        ConductiveTemperatureTransport,
                                        IceWaterThermalEquilibrium
import ClimaSeaIce.SeaIceThermodynamics: initialize_column_interfaces!
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
# `model.phase_transitions` with appropriate parameters. The microscopic
# pure-ice heat capacity enters the Stefan correction to the latent heat:

heat_capacity = 2100 # J kg⁻¹ K⁻¹
phase_transitions = PhaseTransitions(; heat_capacity)

# We set the top ice temperature:

top_temperature = -10 # °C
top_heat_boundary_condition = PrescribedTemperature(top_temperature)

# ## Constructing the thermodynamics
#
# We construct the ice thermodynamics using a simple slab sea ice representation:

ice_thermodynamics = SlabThermodynamics(grid;
                                        internal_heat_flux,
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

model = SeaIceModel(grid; ice_thermodynamics, phase_transitions, sea_ice_density=900, bottom_heat_flux)

# Note that the default bottom heat boundary condition for `SlabThermodynamics`
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

# ## The same bucket with the resolved column thermodynamics
#
# `ColumnEnergyThermodynamics` reads the same `model.external_heat_fluxes` as the slab, so we can freeze an
# identical bucket with a vertically-resolved column and compare. We give the column a vertical grid (a
# `SeaIceColumnDiscretization`) and the same cold lid, internal conductivity, and frazil flux. The column needs a
# tiny non-zero seed thickness for its moving grid, but it still nucleates from open water (ℵ = 0) just like the
# slab. The two freeze nearly identically.

column_grid = RectilinearGrid(size=16, z=SeaIceColumnDiscretization((0, 1)), topology=(Flat, Flat, Bounded))

column_relation = QuadraticLiquidusEnergyRelation(eltype(column_grid); phase_transitions)

column_thermodynamics = prescribed_salinity_enthalpy_thermodynamics(column_grid;
    relation = column_relation,
    salinity_profile = 0,
    energy_transport = ConductiveTemperatureTransport(; conductivity),
    heat_boundary_conditions = (top = top_heat_boundary_condition,
                                bottom = IceWaterThermalEquilibrium(salinity = 0)))

column_model = SeaIceModel(column_grid;
    ice_thermodynamics = column_thermodynamics,
    phase_transitions, sea_ice_density = 900,
    top_heat_flux = 0, bottom_heat_flux)

set!(column_model, h=0.01, ℵ=0)
initialize_column_interfaces!(column_grid, column_model.ice_thickness)
set!(column_thermodynamics; bulk_salinity = 0, temperature = -5)

column_simulation = Simulation(column_model, Δt=10minute, stop_time=10days)

column_timeseries = []
accumulate_column(sim) = push!(column_timeseries, (time(sim),
                                                   first(sim.model.ice_thickness),
                                                   first(sim.model.ice_concentration)))
column_simulation.callbacks[:save] = Callback(accumulate_column)
run!(column_simulation)

tc = [datum[1] for datum in column_timeseries]
hc = [datum[2] for datum in column_timeseries]
ℵc = [datum[3] for datum in column_timeseries]

# All that's left, really, is to put those `lines!` in an `Axis`:

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(size=(1600, 700))

axh = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Ice thickness (cm)")
axℵ = Axis(fig[1, 2], xlabel="Time (days)", ylabel="Ice concentration (-)")
axV = Axis(fig[1, 3], xlabel="Ice volume (cm)", ylabel="Freezing rate (μm s⁻¹)")

lines!(axh, t  ./ day, 1e2 .* h,  label="slab")
lines!(axh, tc ./ day, 1e2 .* hc, linestyle=:dash, label="resolved column")
lines!(axℵ, t  ./ day, ℵ,  label="slab")
lines!(axℵ, tc ./ day, ℵc, linestyle=:dash, label="resolved column")
lines!(axV, 1e2 .* V[1:end-1], 1e6 .* dVdt)
axislegend(axh, position=:lt)

save("freezing_bucket.png", fig)
nothing # hide

# ![](freezing_bucket.png)
#
# If you want more ice, you can increase `simulation.stop_time` and
# `run!(simulation)` again (or just re-run the whole script).
