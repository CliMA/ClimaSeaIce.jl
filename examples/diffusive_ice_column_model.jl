using ClimaSeaIce: ThermodynamicIceModel, MolecularDiffusivity
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators: Δzᶜᶜᶠ
using Oceananigans.Grids: znode
using GLMakie

# Build a grid with 10 cm resolution
grid = RectilinearGrid(size=30, z=(-1, 0), topology=(Flat, Flat, Bounded))

# Set up a simple problem and build the ice model
atmosphere_temperature      = -10  # ᵒC
ocean_temperature           = 1    # ᵒC
closure = MolecularDiffusivity(grid)
model = ThermodynamicIceModel(; grid, closure, atmosphere_temperature, ocean_temperature)

# Initialize and run
set!(model, T=ocean_temperature)

H = model.state.H
@show H

# We're using explicit time stepping. The CFL condition is
#
#   Δt < Δz² / κ ≈ 0.1² / 1e-6 ≈ 1e4,

κ = 1e-5
Δz = Δzᶜᶜᶠ(1, 1, 1, grid)
Δt = 0.1 * Δz^2 / κ
simulation = Simulation(model; Δt)

const c = Center()

ht = Float64[]
th = Float64[]

function compute_ice_thickness(sim; melting_temperature=0.0)
    Nz = size(grid, 3)
    T = sim.model.state.T # model temperature
    Ti = interior(T, 1, 1, :)
    kw = searchsortedlast(Ti, melting_temperature, rev=true)
    ki = kw + 1 # ice index

    # Linearly interpolate to obtain thickness
    Δz = Δzᶜᶜᶠ(1, 1, ki, grid)
    zw = znode(1, 1, kw, grid, c, c, c)
    dTdz = (T[1, 1, ki] - T[1, 1, kw]) / Δz
    z = zw + (melting_temperature - T[1, 1, kw]) / dTdz
    h = -z # thickness

    # Store result
    push!(ht, h)
    push!(th, time(sim))

    return nothing
end

simulation.callbacks[:thickness] = Callback(compute_ice_thickness, IterationInterval(1))

tt = [] 
Tt = [] 
Ht = []
ϕt = []
κt = []

function grab_profiles!(sim)
    # Fields of interest
    T = sim.model.state.T
    H = sim.model.state.H
    ϕ = sim.model.state.ϕ
    κ = sim.model.closure.κ

    # Extract interior data (excluding halos)
    Ti = interior(T, 1, 1, :)  
    Hi = interior(H, 1, 1, :)
    ϕi = interior(ϕ, 1, 1, :)
    κi = interior(κ, 1, 1, :)

    # Copy data into time-series vectors
    push!(tt, time(sim))
    push!(Tt, deepcopy(Ti))
    push!(Ht, deepcopy(Hi))
    push!(ϕt, deepcopy(ϕi))
    push!(κt, deepcopy(κi))

    return nothing
end

simulation.callbacks[:grabber] = Callback(grab_profiles!, SpecifiedTimes(10minutes, 30minutes, 1hour))
simulation.stop_time = 1hour
run!(simulation)

# Make a plot

fig = Figure()

axT = Axis(fig[1, 1], xlabel="Temperature (ᵒC)", ylabel="z (m)")
axϕ = Axis(fig[1, 2], xlabel="Solid fraction", ylabel="z (m)")
axκ = Axis(fig[1, 3], xlabel="Diffusivity", ylabel="z (m)")
axh = Axis(fig[2, 1:3], xlabel="Time (hours)", ylabel="Ice thickness (m)")

z = znodes(model.state.T)

for n = 1:length(tt)
    tn = tt[n]
    Tn = Tt[n]
    ϕn = ϕt[n]
    κn = κt[n]
    label = "t = " * prettytime(tn)
    lines!(axT, Tn, z; label)
    lines!(axϕ, ϕn, z; label)
    lines!(axκ, κn, z; label)
end

lines!(axh, th ./ hour, ht)

axislegend(axT, position=:lb)
display(fig)

