# # Sea ice advected by an atmospheric anticyclone
#
#
#
#
using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using Printf
using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies

# The experiment found in the paper: 
# Simulating Linear Kinematic Features in Viscous-Plastic Sea Ice Models 
# on Quadrilateral and Triangular Grids With Different Variable Staggering

arch = CPU()

L  = 512kilometers
𝓋ₒ = 0.01 # m / s maximum ocean speed
𝓋ₐ = 30.0 # m / s maximum atmospheric speed modifier
Cᴰ = 1.2e-3 # Atmosphere - sea ice drag coefficient
ρₐ = 1.3  # kg/m³
ρₒ = 1025 # kg/m³

# 2 km domain
grid = RectilinearGrid(arch;
                       size = (256, 256), 
                          x = (0, L), 
                          y = (0, L), 
                   topology = (Bounded, Bounded, Flat))

#####
##### Setup atmospheric and oceanic forcing
#####

# Constant ocean velocities corresponding to a cyclonic eddy
Uₒ = XFaceField(grid)
Vₒ = YFaceField(grid)

set!(Uₒ, (x, y) -> 𝓋ₒ * (2y - L) / L)
set!(Vₒ, (x, y) -> 𝓋ₒ * (L - 2x) / L)

struct ExplicitOceanSeaIceStress{U, V, C}
    u    :: U
    v    :: V
    ρₒCᴰ :: C
end

import ClimaSeaIce.SeaIceMomentumEquations: τx, τy

@inline function τx(i, j, k, grid, τ::ExplicitOceanSeaIceStress, clock, fields) 
    Δu = @inbounds fields.u[i, j, k] - τ.u[i, j, k]
    Δv = ℑxyᶠᶜᵃ(i, j, k, grid, τ.v) - ℑxyᶠᶜᵃ(i, j, k, grid, fields.v) 
    return τ.ρₒCᴰ * sqrt(Δu^2 + Δv^2) * Δu
end

@inline function τy(i, j, k, grid, τ::ExplicitOceanSeaIceStress, clock, fields) 
    Δv = @inbounds fields.v[i, j, k] - τ.v[i, j, k] 
    Δu = ℑxyᶜᶠᵃ(i, j, k, grid, τ.u) - ℑxyᶜᶠᵃ(i, j, k, grid, fields.u) 
    return τ.ρₒCᴰ * sqrt(Δu^2 + Δv^2) * Δv
end

τₒ = ExplicitOceanSeaIceStress(U₀, V₀, 5.5)

# Atmosphere - sea ice stress
τᵤₐ = Field(ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Uₐ)
τᵥₐ = Field(ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Vₐ)

# Atmospheric velocities corresponding to an anticyclonic eddy moving north-east
@inline center(t) = 256kilometers + 51.2kilometers * t / 86400
@inline radius(x, y, t)  = sqrt((x - center(t))^2 + (y - center(t))^2)
@inline speed(x, y, t)   = 1 / 100 * exp(- radius(x, y, t) / 100kilometers)

@inline ua_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (cosd(72) * (x - center(t)) + sind(72) * (y - center(t))) / 1000
@inline va_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (cosd(72) * (y - center(t)) - sind(72) * (x - center(t))) / 1000

# Initialize the stress at time t = 0
set!(Uₐ, (x, y) -> ua_time(x, y, 0))
set!(Vₐ, (x, y) -> va_time(x, y, 0))
compute!(τᵤₐ)
compute!(τᵥₐ)

#####
##### Numerical details
#####

# We use an elasto-visco-plastic rheology and WENO seventh order 
# for advection of h and ℵ
momentum_equations = SeaIceMomentumEquation(grid; 
                                            coriolis = FPlane(f=1e-4),
                                            rheology = ElastoViscoPlasticRheology(),
                                            solver   = SplitExplicitSolver(substeps=120))
advection = WENO(; order = 7)

u_bcs = FieldBoundaryConditions(north = ValueBoundaryCondition(0),
                                south = ValueBoundaryCondition(0))

v_bcs = FieldBoundaryConditions(west = ValueBoundaryCondition(0),
                                east = ValueBoundaryCondition(0))

# Define the model!
model = SeaIceModel(grid; 
                    top_momentum_stress = (u = τᵤₐ, v = τᵥₐ),
                    bottom_momentum_stress = (u = τₒ, v = τₒ),
                    ice_dynamics = momentum_equations,
                    ice_thermodynamics = nothing, # No thermodynamics here
                    advection,
                    timestepper = :QuasiAdamsBashforth2,
                    boundary_conditions = (u = u_bcs, v = v_bcs))

# Initial height field with perturbations around 0.3 m
h₀(x, y) = 0.3 + 0.005 * (sin(60 * x / 1000kilometers) + sin(30 * y / 1000kilometers))

# We start with a concentration of ℵ = 1
set!(model, h = h₀)
set!(model, ℵ = 1)

#####
##### Setup the simulation
#####

# run the model for 2 days
simulation = Simulation(model, Δt = 2minutes, stop_time = 2days)

# Remember to evolve the wind stress field in time!
function compute_wind_stress(sim)
    time = sim.model.clock.time
    @inline ua(x, y) = ua_time(x, y, time)
    @inline va(x, y) = va_time(x, y, time)
    set!(Uₐ, ua)
    set!(Vₐ, va)

    compute!(τᵤₐ)
    compute!(τᵥₐ)

    return nothing
end

simulation.callbacks[:top_stress] = Callback(compute_wind_stress, IterationInterval(1))

# Container to hold the data
htimeseries = []
ℵtimeseries = []
utimeseries = []
vtimeseries = []

## Callback function to collect the data from the `sim`ulation
function accumulate_timeseries(sim)
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    u = sim.model.velocities.u
    v = sim.model.velocities.v
    push!(htimeseries, deepcopy(interior(h)))
    push!(ℵtimeseries, deepcopy(interior(ℵ)))
    push!(utimeseries, deepcopy(interior(u)))
    push!(vtimeseries, deepcopy(interior(v)))
end

wall_time = [time_ns()]

function progress(sim) 
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    u = sim.model.velocities.u
    v = sim.model.velocities.v

    hmax = maximum(interior(h))
    ℵmin = minimum(interior(ℵ))
    umax = maximum(interior(u)), maximum(interior(v))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e), max(h): %.2f, min(ℵ): %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., hmax, ℵmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(5))
simulation.callbacks[:save]     = Callback(accumulate_timeseries, IterationInterval(5))

run!(simulation)

using CairoMakie

# Visualize!
Nt = length(htimeseries)
iter = Observable(1)

hi = @lift(htimeseries[$iter][:, :, 1])
ℵi = @lift(ℵtimeseries[$iter][:, :, 1])
ui = @lift(utimeseries[$iter][:, :, 1])
vi = @lift(vtimeseries[$iter][:, :, 1])

fig = Figure()
ax = Axis(fig[1, 1], title = "sea ice thickness")
heatmap!(ax, hi, colormap = :magma,         colorrange = (0.23, 0.37))

ax = Axis(fig[1, 2], title = "sea ice concentration")
heatmap!(ax, ℵi, colormap = Reverse(:deep), colorrange = (0.75, 1))

ax = Axis(fig[2, 1], title = "zonal velocity")
heatmap!(ax, ui, colorrange = (-0.1, 0.1))

ax = Axis(fig[2, 2], title = "meridional velocity")
heatmap!(ax, vi, colorrange = (-0.1, 0.1))

CairoMakie.record(fig, "sea_ice_dynamics.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
    @info "doing iter $i"
end
