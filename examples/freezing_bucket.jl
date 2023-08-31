using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using GLMakie

# Create a 0D grid
grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

model = SlabSeaIceModel(grid;
                        internal_thermal_flux = ConductiveFlux(conductivity=2.0),
                        surface_temperature = - 10, 
                        top_thermal_boundary_condition = PrescribedTemperature())

simulation = Simulation(model, Δt=10minute, stop_time=100days)
set!(model.ice_thickness, 0.01)

# Accumulate data
timeseries = []

function accumulate_timeseries(sim)
    h = model.ice_thickness
    push!(timeseries, (time(sim), first(h)))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

run!(simulation)

# Extract and visualize data
t = map(ts -> ts[1], timeseries)
h = map(ts -> ts[2], timeseries)

# Compute the velocity of the interface
Δt = simulation.Δt
dhdt = @. (h[2:end] - h[1:end-1]) / Δt

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(1200, 600))

axh = Axis(fig[1, 1], xlabel="Time (hours)", ylabel="Ice thickness (m)")
axd = Axis(fig[2, 1], xlabel="Ice thickness (m)", ylabel="Freezing rate (m s⁻¹)")

lines!(axh, t ./ hour, h)
lines!(axd, h[1:end-1], dhdt)

display(fig)

