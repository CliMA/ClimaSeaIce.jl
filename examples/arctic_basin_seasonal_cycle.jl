# # Arctic basin seasonal cycle example
#
# This example reproduces results from Semtner (1976), which simulates the seasonal
# cycle of sea ice in the Arctic basin. The model uses climatological forcing data
# from Fletcher (1965) for shortwave radiation, longwave radiation, sensible heat
# flux, and latent heat flux. This example demonstrates
#
#   * How to use time-varying climatological forcing data.
#   * How to set up seasonal cycles with `FieldTimeSeries` and cyclical indexing.
#   * How to combine multiple heat flux components.
#   * How to simulate multi-year seasonal cycles.
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
# We use a single grid point to model a horizontally uniform ice column:

using Oceananigans
using Oceananigans.Units
using Oceananigans.Units: Time
using ClimaSeaIce
using CairoMakie

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# ## Climatological forcing data
#
# The forcing data comes from Semtner (1976, table 1), originally tabulated by
# Fletcher (1965). Note: these are in kcal, which was deprecated in the ninth
# General Conference on Weights and Measures in 1948. We convert these to Joules
# (and then to Watts) below.
#
# Month:        Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep   Oct    Nov    Dec

tabulated_shortwave = - [   0,      0,    1.9,    9.9,   17.7,   19.2,   13.6,    9.0,    3.7,   0.4,      0,      0] .* 1e4 # kcal m⁻²
tabulated_longwave  = - [10.4,   10.3,   10.3,   11.6,   15.1,   18.0,   19.1,   18.7,   16.5,  13.9,   11.2,   10.9] .* 1e4 # kcal m⁻²
tabulated_sensible  = - [1.18,   0.76,   0.72,   0.29,  -0.45,  -0.39,  -0.30,  -0.40,  -0.17,   0.1,   0.56,   0.79] .* 1e4 # kcal m⁻²
tabulated_latent    = - [   0,  -0.02,  -0.03,   -0.09, -0.46,  -0.70,  -0.64,  -0.66,  -0.39,  -0.19, -0.01,  -0.01] .* 1e4 # kcal m⁻²

# We assume every month has 30 days and set up the time array:

Nmonths = 12
month_days = 30
year_days = month_days * Nmonths
times_days = 15:30:(year_days - 15)
times = times_days .* day # times in seconds

# Sea ice emissivity:

ϵ = 1

# Convert fluxes to the right units (W m⁻²):

kcal_to_joules = 4184
tabulated_shortwave .*= kcal_to_joules / (month_days * days) 
tabulated_longwave  .*= kcal_to_joules / (month_days * days) .* ϵ
tabulated_sensible  .*= kcal_to_joules / (month_days * days)
tabulated_latent    .*= kcal_to_joules / (month_days * days)

# Let's visualize the forcing data:

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Heat flux (W m⁻²)")
lines!(ax, times ./ day, tabulated_shortwave, label="Shortwave")
lines!(ax, times ./ day, tabulated_longwave, label="Longwave")
lines!(ax, times ./ day, tabulated_sensible, label="Sensible")
lines!(ax, times ./ day, tabulated_latent, label="Latent")
axislegend(ax)
current_figure() #hide

# ## Creating time-varying flux functions
#
# We create `FieldTimeSeries` objects with cyclical indexing so the forcing
# repeats annually:

Rs = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Cyclical())
Rl = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Cyclical())
Qs = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Cyclical())
Ql = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Cyclical())

for (i, time) in enumerate(times)
    set!(Rs[i], tabulated_shortwave[i:i])
    set!(Rl[i], tabulated_longwave[i:i])
    set!(Qs[i], tabulated_sensible[i:i])
    set!(Ql[i], tabulated_latent[i:i])
end

# We define flux functions that linearly interpolate the time series:

@inline function linearly_interpolate_flux(i, j, grid, Tₛ, clock, model_fields, flux)
    t  = Time(clock.time)
    return flux[i, j, 1, t]
end

# For shortwave radiation, we account for surface albedo:

@inline function linearly_interpolate_solar_flux(i, j, grid, Tₛ, clock, model_fields, flux)
    Q = linearly_interpolate_flux(i, j, grid, Tₛ, clock, model_fields, flux)
    α = ifelse(Tₛ < -0.1, 0.75, 0.64) # albedo depends on surface temperature
    return Q * (1 - α)
end

# NOTE: Semtner (1976) uses a wrong value for the Stefan-Boltzmann constant
# (roughly 2% higher) to match the results of Maykut & Untersteiner (1965).
# We use the same value here for comparison purposes:

σ = 5.67e-8 * 1.02 # Wrong value, but used for comparison!

# We create flux functions for each component:

Q_shortwave = FluxFunction(linearly_interpolate_solar_flux, parameters=Rs)
Q_longwave  = FluxFunction(linearly_interpolate_flux,       parameters=Rl)
Q_sensible  = FluxFunction(linearly_interpolate_flux,       parameters=Qs)
Q_latent    = FluxFunction(linearly_interpolate_flux,       parameters=Ql)
Q_emission  = RadiativeEmission(emissivity=ϵ, stefan_boltzmann_constant=σ)

# We combine all top heat fluxes:

top_heat_flux = (Q_shortwave, Q_longwave, Q_sensible, Q_latent, Q_emission)

# ## Building the model
#
# We assemble the model and initialize it with 30 cm of ice at full concentration:

model = SeaIceModel(grid; top_heat_flux)
set!(model, h=0.3, ℵ=1) # Start from 30 cm of ice and full concentration

# ## Running the simulation
#
# We run the simulation for 30 years with an 8-hour time step:

simulation = Simulation(model, Δt=8hours, stop_time=30 * 360days)

# ## Collecting data
#
# We set up a callback to accumulate time series data:

series = []

function accumulate_timeseries(sim)
    T = model.ice_thermodynamics.top_surface_temperature
    h = model.ice_thickness
    ℵ = model.ice_concentration
    Qe = model.external_heat_fluxes.top
    Qe = ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions.getflux(Qe, 1, 1, grid, first(T), model.clock, model.ice_thickness)
    push!(series, (time(sim), first(h), first(T), first(ℵ), Qe))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

run!(simulation)

# ## Visualizing the results
#
# We extract the time series data and visualize the evolution:

t = [datum[1] for datum in series]
h = [datum[2] for datum in series]
T = [datum[3] for datum in series]
ℵ = [datum[4] for datum in series]
Q = [datum[5] for datum in series]

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(size=(1000, 1200))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top temperature (ᵒC)")
axh = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")
axℵ = Axis(fig[3, 1], xlabel="Time (days)", ylabel="Ice concentration (-)")
axQ = Axis(fig[4, 1], xlabel="Time (days)", ylabel="Top heat flux (W m⁻²)")

lines!(axT, t ./ day, T)
lines!(axh, t ./ day, h)
lines!(axℵ, t ./ day, ℵ)
lines!(axQ, t ./ day, Q)

current_figure() #hide

save("ice_timeseries.png", fig)
nothing # hide

# ![](ice_timeseries.png)
#
# The results show the seasonal cycle of sea ice, with ice growing in winter and
# melting in summer. After several years, the model reaches a quasi-equilibrium
# seasonal cycle.
