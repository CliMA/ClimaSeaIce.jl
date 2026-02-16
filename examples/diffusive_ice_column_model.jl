# # Diffusive ice column model example
#
# This example demonstrates the use of `EnthalpyMethodSeaIceModel`, which resolves
# the vertical structure of sea ice using an enthalpy-based formulation. The model
# includes internal heat conduction and tracks the phase transitions between ice and
# water within the ice column. This example demonstrates how to:
#
#   * set up an `EnthalpyMethodSeaIceModel` with molecular diffusivity,
#   * prescribe time-varying boundary conditions,
#   * compute ice thickness from the temperature profile,
#   * visualize the evolution of temperature, enthalpy, porosity, and diffusivity.
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
# We build a one-dimensional vertical grid with 20 grid points and 10 cm resolution:

using ClimaSeaIce.EnthalpyMethodSeaIceModels: EnthalpyMethodSeaIceModel, MolecularDiffusivity
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators: Δzᶜᶜᶠ
using Oceananigans.Grids: znode
using CairoMakie

grid = RectilinearGrid(size=20, z=(-1, 0), topology=(Flat, Flat, Bounded))

# ## Model configuration
#
# We set up a simple problem with molecular diffusivity. The diffusivity differs
# for ice and water phases:

closure = MolecularDiffusivity(grid, κ_ice=1e-5, κ_water=1e-6)

# ## Temperature boundary conditions
#
# We define time-varying boundary conditions for the top (air-ice interface) and
# bottom (ice-ocean interface) of the ice column:

initial_air_ice_temperature = -5 # °C
top_T_amplitude = 5 # °C
top_T_slope = -0.5 / day # °C s⁻¹

# Information about ocean cooling:

initial_ice_ocean_temperature = 1.1 # °C
bottom_T_slope = -0.1 / day # °C s⁻¹

# The top temperature oscillates with a daily cycle and has a long-term cooling trend:

air_ice_temperature(x, y, t) = top_T_slope * t + top_T_amplitude * sin(2π*t/day) + initial_air_ice_temperature

# The bottom temperature cools linearly:

ice_ocean_temperature(x, y, t) = bottom_T_slope * t + initial_ice_ocean_temperature

# Let's visualize the boundary condition functions:

dt = 10minutes
tf = 10days
t = 0:dt:tf

set_theme!(Theme(fontsize=24, linewidth=3))
fig = Figure()
ax = Axis(fig[1, 1], title="Boundary Conditions", xlabel="Time (s)", ylabel="Temperature (ᵒC)")
lines!(ax, t, air_ice_temperature.(0, 0, t), label="Air-ice surface temperature")
lines!(ax, t, ice_ocean_temperature.(0, 0, t), label="Ice-ocean temperature")
axislegend(ax)

current_figure() #hide

# We create the boundary conditions:

top_T_bc = ValueBoundaryCondition(air_ice_temperature)
bottom_T_bc = ValueBoundaryCondition(ice_ocean_temperature)
T_bcs = FieldBoundaryConditions(top=top_T_bc, bottom=bottom_T_bc)

# ## Building the model
#
# We assemble the model with the grid, closure, and boundary conditions:

model = EnthalpyMethodSeaIceModel(; grid, closure, boundary_conditions=(; T=T_bcs))

# We initialize the temperature field:

set!(model, T=initial_ice_ocean_temperature)

# ## Time stepping
#
# We're using explicit time stepping. The CFL condition is
#
# ```math
# Δt < Δz² / κ ≈ 0.1² / 1e-6 ≈ 1e4
# ```

κ = 1e-5
Δz = Δzᶜᶜᶠ(1, 1, 1, grid)
Δt = 0.1 * Δz^2 / κ
simulation = Simulation(model; Δt)

# ## Diagnostics
#
# We set up diagnostics to compute ice thickness and collect profiles:

const c = Center()

ht = Float64[]
th = Float64[]

function compute_ice_thickness(sim; melting_temperature=-0.1)
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

# We also collect profiles of temperature, enthalpy, porosity, and diffusivity:

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

simulation.callbacks[:grabber] = Callback(grab_profiles!, TimeInterval(1hour))

# ## Running the simulation
#
# We run the simulation for 10 days:

simulation.stop_time = 10days
run!(simulation)

# ## Visualizing the results
#
# We create an interactive visualization with sliders to explore the evolution:

fig = Figure()

axT = Axis(fig[1, 1], xlabel="Temperature (ᵒC)", ylabel="z (m)")
axH = Axis(fig[1, 2], xlabel="Enthalpy (J m⁻³)", ylabel="z (m)")
axϕ = Axis(fig[1, 3], xlabel="Porosity", ylabel="z (m)")
axκ = Axis(fig[1, 4], xlabel="Diffusivity (m² s⁻¹)", ylabel="z (m)")
axh = Axis(fig[2, 1:4], xlabel="Time (hours)", ylabel="Ice thickness (m)")
axq = Axis(fig[3, 1:4], xlabel="Time (hours)", ylabel="Surface temperatures (ᵒC)")

xlims!(axT, -15, 1)
xlims!(axH, -30, 10)
xlims!(axϕ, -0.05, 1.05)

Nt = length(tt)
slider = Slider(fig[4, 1:4], range=1:Nt, startvalue=1)
n = slider.value

z = znodes(model.state.T)

tn = @lift tt[$n]
tnh = @lift tt[$n] / hour
Tn = @lift Tt[$n]
Hn = @lift Ht[$n]
ϕn = @lift ϕt[$n]
κn = @lift κt[$n]
label = @lift "t = " * prettytime(tt[$n])
scatterlines!(axT, Tn, z; label)
scatterlines!(axH, Hn, z; label)
scatterlines!(axϕ, ϕn, z; label)
scatterlines!(axκ, κn, z; label)

lines!(axh, th ./ hour, ht)

hn = @lift begin
    m = searchsortedfirst(th, tt[$n])
    ht[m]
end

scatter!(axh, tnh, hn, markersize=30, color=(:red, 0.5))

lines!(axq, th ./ hour, air_ice_temperature.(0, 0, th), label="Air-ice surface temperature")
lines!(axq, th ./ hour, ice_ocean_temperature.(0, 0, th), label="Ice-ocean temperature")
vlines!(axq, tnh)
axislegend(axq, position=:lb)
axislegend(axT, position=:lb)

current_figure() #hide

# The visualization shows the evolution of the temperature profile, enthalpy, porosity
# (liquid fraction), and thermal diffusivity within the ice column. The ice thickness
# evolves as the column freezes and melts in response to the boundary conditions.
