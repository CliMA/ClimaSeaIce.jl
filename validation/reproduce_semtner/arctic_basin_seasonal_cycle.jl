#####
##### The purpose of this script is to reproduce some results from Semtner (1976)
##### Note that this script probably doesn't work and is a work in progress.
##### Don't hestitate to update this script in a PR!
#####

using Oceananigans
using Oceananigans.Units
using ClimaSeaIce

# Forcings (Semtner 1976, table 1), originally tabulated by Fletcher (1965)
# Note: these are in kcal, which was deprecated in the ninth General Conference on Weights and
# Measures in 1948. We convert these to Joules (and then to Watts) below.
# Month:        Jan    Feb    Mar    Apr    May    Jun    Jul    Aug    Sep   Oct    Nov    Dec
tabulated_shortwave = [   0      0    1.9    9.9   17.7   19.2   13.6    9.0    3.7   0.4      0      0]
tabulated_longwave  = [10.4   10.3   10.3   11.6   15.1   18.0   19.1   18.7   16.5  13.9   11.2   10.9]
tabulated_sensible  = [1.18   0.76   0.72   0.29  -0.45  -0.39  -0.30  -0.40  -0.17   0.1   0.56   0.79]
tabulated_latent    = [   0  -0.02  -0.03  -0.46  -0.70  -0.64  -0.66  -0.39  -0.19 -0.01  -0.01  -3.20]

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

@inline function linearly_interpolate_flux(i, j, grid, Tₛ, clock, model_fields, parameters)
    times = parameters.times
    Q = parameters.flux
    Nt = length(times)
    t = mod(clock.time, 360days)
    n₁, n₂ = index_binary_search(times, t, Nt)

    Q₁ = @inbounds Q[n₁]
    Q₂ = @inbounds Q[n₂]

    ñ = @inbounds (t - times[n₁]) * (n₂ - n₁) / (times[n₂] - times[n₁])

    return Q₁ * (ñ - n₁) + Q₂ * (n₂ - ñ)
end

Q_shortwave = FluxFunction(linearly_interpolate_flux, parameters=(; times, flux=tabulated_shortwave))
Q_longwave  = FluxFunction(linearly_interpolate_flux, parameters=(; times, flux=tabulated_longwave))
Q_sensible  = FluxFunction(linearly_interpolate_flux, parameters=(; times, flux=tabulated_sensible))
Q_latent    = FluxFunction(linearly_interpolate_flux, parameters=(; times, flux=tabulated_latent))
Q_emission  = RadiativeEmission()

# Generate a 0D grid for a single column slab model 
grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

top_heat_flux = (Q_shortwave, Q_longwave, Q_sensible, Q_latent, Q_emission)
model = SlabSeaIceModel(grid; top_heat_flux)
set!(model, h=1)

simulation = Simulation(model, Δt=8hours, stop_time=4 * 360days)

# Accumulate data
timeseries = []

function accumulate_timeseries(sim)
    T = model.top_temperature
    h = model.ice_thickness
    push!(timeseries, (time(sim), first(h), first(T)))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

run!(simulation)

# Extract and visualize data

t = [datum[1] for datum in timeseries]
h = [datum[2] for datum in timeseries]
T = [datum[3] for datum in timeseries]

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(1000, 800))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top temperature (ᵒC)")
axh = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")

lines!(axT, t / day, T)
lines!(axh, t / day, h)

display(fig)

