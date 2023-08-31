using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.ThermalBoundaryConditions: RadiativeEmission
using GLMakie

# Generate a 0D grid for a single column slab model 
grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# Build a model of an ice slab that has internal conductive fluxes
# and that emits radiation from its surface.
model = SlabSeaIceModel(grid; top_external_thermal_fluxes = RadiativeEmission())
set!(model.ice_thickness, 0.01)

simulation = Simulation(model, Δt=1hour, stop_time=40days)

# Accumulate data
timeseries = []

function accumulate_timeseries(sim)
    T = model.surface_temperature
    h = model.ice_thickness
    push!(timeseries, (time(sim), first(h), first(T)))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

run!(simulation)

# Extract and visualize data

t = map(ts -> ts[1], timeseries)
h = map(ts -> ts[2], timeseries)
T = map(ts -> ts[3], timeseries)

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(1000, 800))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top surface temperature (ᵒC)")
axh = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")

lines!(axT, t / day, T)
lines!(axh, t / day, h)

display(fig)

