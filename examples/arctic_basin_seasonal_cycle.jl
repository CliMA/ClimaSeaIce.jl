#####
##### The purpose of this script is to reproduce some results from Semtner (1976)
#####

using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: Time
using ClimaSeaIce

# Generate a 0D grid for a single column slab model 
grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# Forcings (Semtner 1976, table 1), originally tabulated by Fletcher (1965)
# Note: these are in kcal, which was deprecated in the ninth General Conference on Weights and
# Measures in 1948. We convert these to Joules (and then to Watts) below.
# Month:        Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep   Oct    Nov    Dec
tabulated_shortwave = - [   0,      0,    1.9,    9.9,   17.7,   19.2,   13.6,    9.0,    3.7,   0.4,      0,      0] .* 1e4 # to convert to kcal / m²
tabulated_longwave  = - [10.4,   10.3,   10.3,   11.6,   15.1,   18.0,   19.1,   18.7,   16.5,  13.9,   11.2,   10.9] .* 1e4 # to convert to kcal / m²
tabulated_sensible  = - [1.18,   0.76,   0.72,   0.29,  -0.45,  -0.39,  -0.30,  -0.40,  -0.17,   0.1,   0.56,   0.79] .* 1e4 # to convert to kcal / m²
tabulated_latent    = - [   0,  -0.02,  -0.03,   -0.09, -0.46,  -0.70,  -0.64,  -0.66,  -0.39,  -0.19, -0.01,  -0.01] .* 1e4 # to convert to kcal / m²

# Pretend every month is just 30 days
Nmonths = 12
month_days = 30
year_days = month_days * Nmonths
times_days = 15:30:(year_days - 15)
times = times_days .* day # times in seconds

# Convert fluxes to the right units
kcal_to_joules = 4184
tabulated_shortwave .*= kcal_to_joules / (month_days * days) 
tabulated_longwave  .*= kcal_to_joules / (month_days * days)
tabulated_sensible  .*= kcal_to_joules / (month_days * days)
tabulated_latent    .*= kcal_to_joules / (month_days * days)

using GLMakie

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, times, tabulated_shortwave)
lines!(ax, times, tabulated_longwave)
lines!(ax, times, tabulated_sensible)
lines!(ax, times, tabulated_latent)

display(fig)

# Make them into a FieldTimeSeries for better manipulation

Rs = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Cyclical)
Rl = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Cyclical)
Qs = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Cyclical)
Ql = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Cyclical)

for (i, time) in enumerate(times)
    set!(Rs, tabulated_shortwave[i])
    set!(Rl, tabulated_longwave[i])
    set!(Qs, tabulated_sensible[i])
    set!(Ql, tabulated_latent[i])
end

@inline function linearly_interpolate_flux(i, j, grid, Tₛ, clock, model_fields, flux)
    t  = Time(clock.time)
    return flux[i, j, 1, t]
end

@inline function linearly_interpolate_solar_flux(i, j, grid, Tₛ, clock, model_fields, flux)
    Q = linearly_interpolate_flux(i, j, grid, Tₛ, clock, model_fields, flux)
    α = ifelse(Tₛ < -0.1, 0.75, 0.64)
    return Q * (1 - α)
end

Q_shortwave = FluxFunction(linearly_interpolate_solar_flux, parameters=Rs)
Q_longwave  = FluxFunction(linearly_interpolate_flux,       parameters=Rl)
Q_sensible  = FluxFunction(linearly_interpolate_flux,       parameters=Qs)
Q_latent    = FluxFunction(linearly_interpolate_flux,       parameters=Ql)
Q_emission  = RadiativeEmission()

top_heat_flux = (Q_shortwave, Q_longwave, Q_sensible, Q_latent, Q_emission)

model = SeaIceModel(grid; top_heat_flux)
set!(model, h=0.3, ℵ=1) # We start from 300cm of ice and full concentration

simulation = Simulation(model, Δt=1hours, stop_time=2* 360days)

# Accumulate data
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

# Extract and visualize data

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
axQ = Axis(fig[4, 1], xlabel="Time (days)", ylabel="Top heat flux (W/m²)")

lines!(axT, t / day, T)
lines!(axh, t / day, h)
lines!(axℵ, t / day, ℵ)
lines!(axQ, t / day, Q)

display(fig)

