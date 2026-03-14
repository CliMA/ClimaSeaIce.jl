# # Arctic basin seasonal cycle example
#
# This example reproduces results from Semtner (1976), which simulates the seasonal
# cycle of sea ice in the Arctic basin. The model uses climatological forcing data
# from Fletcher (1965) for shortwave radiation, longwave radiation, sensible heat
# flux, and latent heat flux.
#
# We run two side-by-side simulations --- bare ice and snow-covered ice --- to show
# the insulating effect of snow. Snow has a much lower thermal conductivity than ice
# (~0.31 vs ~2 W/m/K), reducing the conductive flux through the slab and therefore
# slowing ice growth.
#
# This example demonstrates how to:
#
#   * use time-varying climatological forcing data,
#   * set up seasonal cycles with `FieldTimeSeries` and cyclical indexing,
#   * combine multiple heat flux components,
#   * add a snow layer with `snow_thermodynamics`,
#   * compare ice growth with and without snow insulation.
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

Rs = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Oceananigans.OutputReaders.Cyclical())
Rl = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Oceananigans.OutputReaders.Cyclical())
Qs = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Oceananigans.OutputReaders.Cyclical())
Ql = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Oceananigans.OutputReaders.Cyclical())

for (i, time) in enumerate(times)
    set!(Rs[i], tabulated_shortwave[i:i])
    set!(Rl[i], tabulated_longwave[i:i])
    set!(Qs[i], tabulated_sensible[i:i])
    set!(Ql[i], tabulated_latent[i:i])
end

# We define flux functions that linearly interpolate the time series:

@inline function linearly_interpolate_flux(i, j, grid, Ts, clock, model_fields, flux)
    t  = Time(clock.time)
    return flux[i, j, 1, t]
end

# For shortwave radiation, we account for surface albedo:

@inline function linearly_interpolate_solar_flux(i, j, grid, Ts, clock, model_fields, flux)
    Q = linearly_interpolate_flux(i, j, grid, Ts, clock, model_fields, flux)
    α = ifelse(Ts < -0.1, 0.75, 0.64) # albedo depends on surface temperature
    return Q * (1 - α)
end

# !!! note "Slightly wrong value for Stefan-Boltzmann constant"
#
#     Semtner (1976) uses a wrong value for the Stefan-Boltzmann constant
#     (roughly 2% higher) to match the results of Maykut & Untersteiner (1965).
#     We use here the same value here for comparison purposes.

σ = 5.67e-8 * 1.02 # Wrong value, but used for comparison!

# We create flux functions for each component:

Q_shortwave = FluxFunction(linearly_interpolate_solar_flux, parameters=Rs)
Q_longwave  = FluxFunction(linearly_interpolate_flux,       parameters=Rl)
Q_sensible  = FluxFunction(linearly_interpolate_flux,       parameters=Qs)
Q_latent    = FluxFunction(linearly_interpolate_flux,       parameters=Ql)
Q_emission  = RadiativeEmission(emissivity=ϵ, stefan_boltzmann_constant=σ)

top_heat_flux = (Q_shortwave, Q_longwave, Q_sensible, Q_latent, Q_emission)

# ## Building the bare-ice model
#
# We assemble the bare-ice model and initialize it with 30 cm of ice at full concentration:

model_bare = SeaIceModel(grid; top_heat_flux)
set!(model_bare, h=0.3, ℵ=1)

# ## Snow thermodynamics
#
# Snow has a much lower thermal conductivity than ice. When it accumulates on top
# of sea ice, it adds thermal resistance in series with the ice: R = hs/ks + hi/ki.
# This reduces the conductive flux that drives ice growth at the bottom.
#
# We set up the snow layer with physical properties for Arctic snow:

snow_thermodynamics = SlabSnowThermodynamics(grid)

# We prescribe a constant snowfall rate of about 3e-6 kg/m²/s, corresponding
# to roughly 30 cm/year of snow accumulation at ρs = 330 kg/m³.
# In summer, the excess surface energy melts snow faster than it accumulates,
# so the snow cover follows a natural seasonal cycle.

snow_precipitation = 3e-6 # kg m⁻² s⁻¹

# ## Building the snow-covered model
#
# When `snow_thermodynamics` is provided, the `SeaIceModel` constructor
# automatically builds the layered coupling: the snow layer sits on top of
# ice, drives the surface temperature solve using the combined conductive
# flux, and computes the interface temperature between snow and ice.

model_snow = SeaIceModel(grid; top_heat_flux, snow_thermodynamics, snow_precipitation)
set!(model_snow, h=0.3, ℵ=1, hs=0.05) # Start with 5 cm of snow

# ## Running the simulations
#
# We run both simulations for 10 years with an 8-hour time step:

Δt = 8hours
stop_time = 10 * 360days

# Bare-ice simulation:

simulation_bare = Simulation(model_bare, Δt=Δt, stop_time=stop_time)

series_bare = []

function accumulate_bare(sim)
    T = sim.model.ice_thermodynamics.top_surface_temperature
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    push!(series_bare, (time(sim), first(h), first(T), first(ℵ)))
end

simulation_bare.callbacks[:save] = Callback(accumulate_bare)
run!(simulation_bare)

# Snow-covered simulation:

simulation_snow = Simulation(model_snow, Δt=Δt, stop_time=stop_time)

series_snow = []

function accumulate_snow(sim)
    Tu = sim.model.snow_thermodynamics.top_surface_temperature
    h  = sim.model.ice_thickness
    ℵ  = sim.model.ice_concentration
    hs = sim.model.snow_thickness
    push!(series_snow, (time(sim), first(h), first(Tu), first(ℵ), first(hs)))
end

simulation_snow.callbacks[:save] = Callback(accumulate_snow)
run!(simulation_snow)

# ## Visualizing the results
#
# We extract the time series data and compare:

t_bare = [d[1] for d in series_bare]
h_bare = [d[2] for d in series_bare]
T_bare = [d[3] for d in series_bare]

t_snow = [d[1] for d in series_snow]
h_snow = [d[2] for d in series_snow]
T_snow = [d[3] for d in series_snow]
hs     = [d[5] for d in series_snow]

set_theme!(Theme(fontsize=24, linewidth=3))

fig = Figure(size=(1200, 1200))

axT = Axis(fig[1, 1], ylabel="Surface temperature (°C)")
lines!(axT, t_bare ./ day, T_bare, label="Bare ice")
lines!(axT, t_snow ./ day, T_snow, label="Snow-covered", linestyle=:dash)
axislegend(axT, position=:rb)

axh = Axis(fig[2, 1], ylabel="Ice thickness (m)")
lines!(axh, t_bare ./ day, h_bare, label="Bare ice")
lines!(axh, t_snow ./ day, h_snow, label="Snow-covered", linestyle=:dash)
axislegend(axh, position=:rb)

axs = Axis(fig[3, 1], ylabel="Snow thickness (m)")
lines!(axs, t_snow ./ day, hs, color=Makie.wong_colors()[2])

axd = Axis(fig[4, 1], xlabel="Time (days)", ylabel="Thickness difference (m)")
lines!(axd, t_bare ./ day, h_bare .- h_snow[1:length(h_bare)], color=:black)
hlines!(axd, [0], color=:gray, linestyle=:dot, linewidth=1)

save("ice_timeseries.png", fig)
nothing # hide

# ![](ice_timeseries.png)
#
# The results show that snow insulation has a significant effect on the seasonal
# cycle. The snow-covered ice:
#
#   * **grows more slowly in winter** because the snow layer reduces the conductive
#     flux from the cold surface to the ice-water interface,
#   * **has a warmer surface temperature** because the same heat flux produces a
#     larger temperature drop across the more insulating snow+ice slab,
#   * **accumulates snow in autumn/winter** and **melts it in spring/summer** as
#     the surface energy balance turns positive.
#
# The bottom panel shows the thickness difference: positive means bare ice is
# thicker, confirming that snow insulation reduces net ice growth.
