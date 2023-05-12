using ClimaSeaIce: ThermodynamicIceModel, ConstantIceDiffusivity
using Oceananigans
using Oceananigans.Units
using GLMakie

# Build a grid with 10 cm resolution
grid = RectilinearGrid(size=10, z=(-1, 0), topology=(Flat, Flat, Bounded))

# Set up a simple problem and build the ice model
atmosphere_temperature      = -10  # ᵒC
ocean_temperature           = 0    # ᵒC
ice_temperature_diffusivity = 1e-6 # m² s⁻¹
closure = ConstantIceDiffusivity(ice_temperature_diffusivity)
model = ThermodynamicIceModel(; grid, closure, atmosphere_temperature, ocean_temperature)

# Initialize and run
model.temperature .= 0 

# We're using explicit time stepping. The CFL condition is
#
#   Δt < Δz² / κ ≈ 0.1² / 1e-6 ≈ 1e4,
#
# So Δt=1minute is well within CFL.
simulation = Simulation(model, Δt=1minute)

fig = Figure()
ax = Axis(fig[1, 1], xlabel="Temperature (ᵒC)", ylabel="z (m)")

for stop_time in (10minutes, 30minutes, 1hour)
    simulation.stop_time = stop_time
    run!(simulation)

    T = model.temperature
    Tn = interior(T, 1, 1, :)
    z = znodes(T)

    label = "t = " * prettytime(model.clock.time)
    lines!(ax, Tn, z; label)
end

axislegend(ax, position=:lb)
display(fig)

