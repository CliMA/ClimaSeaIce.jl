# # Melting in spring example
#
# This example simulates the melting of relatively thick sea ice in spring when
# the sun is shining. The ice is subject to solar insolation and sensible heat
# fluxes from the atmosphere. Different grid cells show how the ice melts at
# different rates depending on the amount of solar insolation they receive.
# This example demonstrates how to:
#
#   * set up a one-dimensional model with multiple grid cells,
#   * prescribe spatially varying solar insolation,
#   * use `FluxFunction` for parameterized heat fluxes,
#   * visualize the evolution of multiple ice columns.
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
# We generate a one-dimensional grid with 4 grid cells to model different ice
# columns subject to different solar insolation:

using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: RadiativeEmission
using CairoMakie

grid = RectilinearGrid(size=4, x=(0, 1), topology=(Periodic, Flat, Flat))

# ## Top boundary conditions
#
# We prescribe different solar insolation values for each grid cell, ranging
# from -600 to -1200 W m⁻² (negative values indicate downward/into the ice):

solar_insolation = [-600, -800, -1000, -1200] # W m⁻²
solar_insolation = reshape(solar_insolation, (4, 1, 1))

# The ice also emits longwave radiation from its top surface:

outgoing_radiation = RadiativeEmission()

# ## Sensible heat flux parameterization
#
# The sensible heat flux from the atmosphere is represented by a `FluxFunction`.
# We define the parameters for the bulk formula:

parameters = (
    transfer_coefficient     = 1e-3,  # unitless
    atmosphere_density       = 1.225, # kg m⁻³
    atmosphere_heat_capacity = 1004,  # J kg⁻¹ K⁻¹
    atmosphere_temperature   = -5,    # °C
    atmosphere_wind_speed    = 5      # m s⁻¹
)

# The flux is positive (cooling by fluxing heat upward away from the upper surface)
# when the atmosphere temperature is less than the surface temperature:

@inline function sensible_heat_flux(i, j, grid, Tᵤ, clock, fields, parameters)
    Cₛ = parameters.transfer_coefficient
    ρₐ = parameters.atmosphere_density
    cₐ = parameters.atmosphere_heat_capacity
    Tₐ = parameters.atmosphere_temperature
    uₐ = parameters.atmosphere_wind_speed
    ℵ  = fields.ℵ[i, j, 1]

    return Cₛ * ρₐ * cₐ * uₐ * (Tᵤ - Tₐ) * ℵ
end

aerodynamic_flux = FluxFunction(sensible_heat_flux; parameters)

# We combine all top heat fluxes into a tuple:

top_heat_flux = (outgoing_radiation, solar_insolation, aerodynamic_flux)

# ## Building the model
#
# We assemble the model with an ice consolidation thickness of 0.05 m:

model = SeaIceModel(grid;
                    ice_consolidation_thickness = 0.05, # m
                    top_heat_flux)

# We initialize all columns with a 1 m thick slab of ice with 100% ice concentration:

set!(model, h=1, ℵ=1)

# ## Running a simulation
#
# We set up a simulation that runs for 30 days with a 10-minute time step:

simulation = Simulation(model, Δt=10minute, stop_time=30days)

# ## Collecting data
#
# We set up a callback to accumulate time series data for all four columns:

timeseries = []

function accumulate_timeseries(sim)
    T = model.ice_thermodynamics.top_surface_temperature
    h = model.ice_thickness
    ℵ = model.ice_concentration
    push!(timeseries, (time(sim),
                       h[1, 1, 1], ℵ[1, 1, 1], T[1, 1, 1],
                       h[2, 1, 1], ℵ[2, 1, 1], T[2, 1, 1],
                       h[3, 1, 1], ℵ[3, 1, 1], T[3, 1, 1],
                       h[4, 1, 1], ℵ[4, 1, 1], T[4, 1, 1]))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

run!(simulation)

# ## Visualizing the results
#
# We extract the time series data for all four columns:

t  = [datum[1]  for datum in timeseries]
h1 = [datum[2]  for datum in timeseries]
ℵ1 = [datum[3]  for datum in timeseries]
T1 = [datum[4]  for datum in timeseries]
h2 = [datum[5]  for datum in timeseries]
ℵ2 = [datum[6]  for datum in timeseries]
T2 = [datum[7]  for datum in timeseries]
h3 = [datum[8]  for datum in timeseries]
ℵ3 = [datum[9]  for datum in timeseries]
T3 = [datum[10] for datum in timeseries]
h4 = [datum[11] for datum in timeseries]
ℵ4 = [datum[12] for datum in timeseries]
T4 = [datum[13] for datum in timeseries]

# And visualize the evolution of ice thickness, concentration, and surface temperature:

set_theme!(Theme(fontsize=18, linewidth=3))

fig = Figure(size=(1000, 800))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top temperature (ᵒC)")
axℵ = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice concentration (-)")
axh = Axis(fig[3, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")

lines!(axT, t ./ day, T1)
lines!(axℵ, t ./ day, ℵ1)
lines!(axh, t ./ day, h1)

lines!(axT, t ./ day, T2)
lines!(axℵ, t ./ day, ℵ2)
lines!(axh, t ./ day, h2)

lines!(axT, t ./ day, T3)
lines!(axℵ, t ./ day, ℵ3)
lines!(axh, t ./ day, h3)

lines!(axT, t ./ day, T4)
lines!(axℵ, t ./ day, ℵ4)
lines!(axh, t ./ day, h4)

save("melting_in_spring.png", fig)
nothing # hide

# ![](melting_in_spring.png)
#
# The results show that ice columns receiving more solar insolation melt faster,
# as expected. The ice concentration and thickness decrease more rapidly for
# columns with higher incoming solar radiation.
