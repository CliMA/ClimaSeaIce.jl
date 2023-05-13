using ClimaSeaIce: ThermodynamicIceModel, ScalarIceDiffusivity
using Oceananigans
using Oceananigans.Units
using GLMakie

# Build a grid with 10 cm resolution
grid = RectilinearGrid(size=10, z=(-1, 0), topology=(Flat, Flat, Bounded))

# Set up a simple problem and build the ice model
atmosphere_temperature      = -10  # ᵒC
ocean_temperature           = 1    # ᵒC
ice_temperature_diffusivity = 1e-6 # m² s⁻¹
closure = ScalarIceDiffusivity(ice_temperature_diffusivity)
model = ThermodynamicIceModel(; grid, closure, atmosphere_temperature, ocean_temperature)

# Initialize and run
set!(model, T=ocean_temperature)

H = model.state.H
@show H

# We're using explicit time stepping. The CFL condition is
#
#   Δt < Δz² / κ ≈ 0.1² / 1e-6 ≈ 1e4,
#
# So Δt=1minute is well within CFL.
simulation = Simulation(model, Δt=1minute)

fig = Figure()
axT = Axis(fig[1, 1], xlabel="Temperature (ᵒC)", ylabel="z (m)")
axϕ = Axis(fig[1, 2], xlabel="Porosity", ylabel="z (m)")

for stop_time in (0, 10minutes, 30minutes, 1hour)
#for stop_time in (0, simulation.Δt)
    simulation.stop_time = stop_time
    run!(simulation)

    T = model.state.T
    H = model.state.H
    ϕ = model.state.ϕ
    Tn = interior(T, 1, 1, :)
    Hn = interior(H, 1, 1, :)
    ϕn = interior(ϕ, 1, 1, :)
    z = znodes(T)

    label = "t = " * prettytime(model.clock.time)
    lines!(axT, Tn, z; label)
    lines!(axϕ, ϕn, z; label)
end

axislegend(ax, position=:lb)
display(fig)

