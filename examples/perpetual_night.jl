using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.ThermalBoundaryConditions: RadiativeEmission
using GLMakie

# Generate a zero-dimensional grid for a single column slab model 
grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# Build a model of an ice slab that has internal conductive fluxes
# and that emits radiation from its top surface.
model = SlabSeaIceModel(grid; top_thermal_flux=RadiativeEmission())
set!(model, h=0.01)

simulation = Simulation(model, Δt=1hour, stop_time=40days)

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

