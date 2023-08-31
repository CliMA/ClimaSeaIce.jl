using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using GLMakie

# Create a single column grid with no grid points in the vertical,
# appropriate for a zero heat capacity model.
Nx = Ny = 1

grid = RectilinearGrid(size = (Nx, Ny),
                       x = (0, 1),
                       y = (0, 1),
                       topology = (Periodic, Periodic, Flat))

internal_thermal_flux = ConductiveFlux(conductivity=1.0)

model = SlabSeaIceModel(grid;
                        surface_temperature = - 10,
                        top_thermal_boundary_condition = PrescribedTemperature())

simulation = Simulation(model, Δt=10minute, stop_time=100days)

model.ice_thickness .= 0.01

# Integrate forward for 1 hour
ht = Float64[]
times = Float64[]

function save_ice_state(sim)
    h = model.ice_thickness
    push!(ht, first(h))
    push!(times, time(sim))
end

simulation.callbacks[:save] = Callback(save_ice_state)

#=
h = model.ice_thickness
simulation.output_writers[:jld2] = JLD2OutputWriter(model, (; h),
                                                    schedule = IterationInterval(1),
                                                    filename = "freezing_bucket.jld2",
                                                    overwrite_existing = true)
=#

run!(simulation)

# Post-processing - calculate dh/dt
Δt = simulation.Δt
dh_dt = (ht[2:end] - ht[1:end-1]) ./ Δt

set_theme!(Theme(fontsize=24, linewidth=4))
fig = Figure(resolution=(1200, 600))
axh = Axis(fig[1, 1], xlabel="Time (hours)", ylabel="Ice thickness (m)")
axd = Axis(fig[1, 2], xlabel="Ice thickness (m)", ylabel="Freezing rate (m s⁻¹)")
lines!(axh, times ./ hour, ht)
lines!(axd, ht[1:end-1], dh_dt)
display(fig)

