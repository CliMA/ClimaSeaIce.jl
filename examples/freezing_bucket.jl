# # Freezing bucket example
#
# A common laboratory experiment freezes an insulated bucket of water from the top
# down, using a metal lid to keep the top of the bucket at some constant, very cold
# temperature. In this example, we simulate such a scenario using the `SeaIceModel`,
# comparing ice growth with and without a snow layer on top.
#
# Snow acts as an insulating blanket: its low thermal conductivity (~0.31 W/m/K vs
# ~2 W/m/K for ice) reduces the conductive heat flux through the slab, slowing ice
# growth. We run two side-by-side simulations to see this effect.
#
# This example demonstrates how to:
#
#   * use `SlabThermodynamics` with prescribed boundary conditions,
#   * configure internal heat conduction,
#   * add a snow layer with `snow_thermodynamics` and `snow_precipitation`,
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
using CairoMakie

# to load ClimaSeaIce functions and objects into our script.
#
# ## An infinitely deep bucket with a single grid point
#
# Perhaps surprisingly, we need just one grid point to model a possibly infinitely
# thick slab of ice with `SeaIceModel`. We would only need more than 1 grid point
# if our boundary conditions vary in the horizontal direction:

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# ## Ice thermodynamics
#
# We build our model of an ice slab freezing into a bucket. We start by defining
# a constant internal `ConductiveFlux` with ice conductivity:

ki = 2 # W m⁻¹ K⁻¹
ice_internal_heat_flux = ConductiveFlux(; conductivity=ki)

# We set the heat capacity and density:

ice_heat_capacity = 2100 # J kg⁻¹ K⁻¹
ice_density = 900 # kg m⁻³
ice_phase_transitions = PhaseTransitions(; heat_capacity=ice_heat_capacity, density=ice_density)

# We set the top ice temperature to a cold -10 °C:

top_temperature = -10 # °C
top_heat_boundary_condition = PrescribedTemperature(top_temperature)

# We construct the ice thermodynamics:

ice_thermodynamics = SlabThermodynamics(grid;
                                        internal_heat_flux = ice_internal_heat_flux,
                                        phase_transitions = ice_phase_transitions,
                                        top_heat_boundary_condition)

# ## Frazil ice formation
#
# We prescribe a frazil ice heat flux that stops when the ice has reached a
# concentration of 1. This heat flux represents the initial ice formation from a
# liquid bucket:

@inline frazil_ice_formation(i, j, grid, Tui, clock, fields) = -(1 - fields.ℵ[i, j, 1]) # W m⁻²

bottom_heat_flux = FluxFunction(frazil_ice_formation)

# ## Model without snow
#
# First, a bare-ice model:

model_bare = SeaIceModel(grid; ice_thermodynamics, bottom_heat_flux)

# ## Snow thermodynamics
#
# Snow has a much lower thermal conductivity than ice, typically around 0.31 W/m/K.
# When snow sits on top of ice, it acts as an insulating layer that reduces the
# total conductive heat flux. The combined thermal resistance is the sum of
# the snow and ice resistances in series: R = hs/ks + hi/ki.
#
# We define the snow layer with its own thermodynamic properties:

snow_thermodynamics = SlabSnowThermodynamics(grid)

# We also add a light snowfall rate (about 3.6 cm/day of snow accumulation):

snow_precipitation = 4e-4 # kg m⁻² s⁻¹

# ## Model with snow
#
# When `snow_thermodynamics` is provided, the `SeaIceModel` constructor
# automatically wires the layered coupling: the snow layer sits on top of
# the ice and drives the surface temperature solve using the combined
# snow+ice conductive flux.

model_snow = SeaIceModel(grid;
    ice_thermodynamics = SlabThermodynamics(grid;
                                            internal_heat_flux = ice_internal_heat_flux,
                                            phase_transitions = ice_phase_transitions,
                                            top_heat_boundary_condition),
    bottom_heat_flux,
    snow_thermodynamics,
    snow_precipitation)

# ## Running both simulations
#
# We freeze the bucket for 30 days with a 10-minute time step:

Δt = 10minute
stop_time = 30days

# First, set up and run the bare-ice simulation:

simulation_bare = Simulation(model_bare, Δt=Δt, stop_time=stop_time)

timeseries_bare = []

function accumulate_bare(sim)
    push!(timeseries_bare, (time(sim),
                            first(sim.model.ice_thickness),
                            first(sim.model.ice_concentration)))
end

simulation_bare.callbacks[:save] = Callback(accumulate_bare)
run!(simulation_bare)

# Now set up and run the snow-covered simulation:

simulation_snow = Simulation(model_snow, Δt=Δt, stop_time=stop_time)

timeseries_snow = []

function accumulate_snow(sim)
    push!(timeseries_snow, (time(sim),
                            first(sim.model.ice_thickness),
                            first(sim.model.ice_concentration),
                            first(sim.model.snow_thickness)))
end

simulation_snow.callbacks[:save] = Callback(accumulate_snow)
run!(simulation_snow)

# ## Visualizing the results
#
# We extract the time series and compare:

t_bare = [d[1] for d in timeseries_bare]
h_bare = [d[2] for d in timeseries_bare]
ℵ_bare = [d[3] for d in timeseries_bare]

t_snow = [d[1] for d in timeseries_snow]
h_snow = [d[2] for d in timeseries_snow]
ℵ_snow = [d[3] for d in timeseries_snow]
hs     = [d[4] for d in timeseries_snow]

# Compute freezing rates (ice volume tendency):

V_bare = h_bare .* ℵ_bare
V_snow = h_snow .* ℵ_snow
dVdt_bare = @. (V_bare[2:end] - V_bare[1:end-1]) / Δt
dVdt_snow = @. (V_snow[2:end] - V_snow[1:end-1]) / Δt

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(size=(1200, 900))

# Ice thickness comparison
axh = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Ice thickness (cm)",
           title="Ice growth: bare vs snow-covered")
lines!(axh, t_bare ./ day, 1e2 .* h_bare, label="Bare ice")
lines!(axh, t_snow ./ day, 1e2 .* h_snow, label="Snow-covered ice", linestyle=:dash)
axislegend(axh, position=:rb)

# Snow thickness
axs = Axis(fig[1, 2], xlabel="Time (days)", ylabel="Snow thickness (cm)",
           title="Snow accumulation")
lines!(axs, t_snow ./ day, 1e2 .* hs, color=:dodgerblue)

# Freezing rate comparison
axf = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Freezing rate (μm s⁻¹)",
           title="Bottom freezing rate")
lines!(axf, t_bare[2:end] ./ day, 1e6 .* dVdt_bare, label="Bare ice")
lines!(axf, t_snow[2:end] ./ day, 1e6 .* dVdt_snow, label="Snow-covered ice", linestyle=:dash)
axislegend(axf, position=:rt)

# Ice concentration
axc = Axis(fig[2, 2], xlabel="Time (days)", ylabel="Ice concentration (-)",
           title="Ice concentration")
lines!(axc, t_bare ./ day, ℵ_bare, label="Bare ice")
lines!(axc, t_snow ./ day, ℵ_snow, label="Snow-covered ice", linestyle=:dash)
axislegend(axc, position=:rb)

save("freezing_bucket.png", fig)
nothing # hide

# ![](freezing_bucket.png)
#
# The snow-covered ice grows significantly slower than bare ice. The snow layer
# acts as an insulating blanket, reducing the conductive heat flux from the cold
# surface to the warm ice-water interface. With less heat escaping through the slab,
# less ice freezes at the bottom.
