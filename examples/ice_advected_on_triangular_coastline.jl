using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using Printf
using CairoMakie
using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies

# A solid block of ice moving against a triangular coastline in a periodic channel

Lx = 512kilometers
Ly = 256kilometers
Nx = 256
Ny = 256

y_max = Ly / 2

arch = GPU()

𝓋ₐ = 10.0   # m / s 
Cᴰ = 1.2e-3 # Atmosphere - sea ice drag coefficient
ρₐ = 1.3    # kg/m³

# 2 km domain
grid = RectilinearGrid(arch; size = (Nx, Ny), 
                                x = (-Lx/2, Lx/2), 
                                y = (0, Ly), 
                             halo = (4, 4),
                         topology = (Periodic, Bounded, Flat))

bottom(x, y) = ifelse(y > y_max, 0, 
               ifelse(abs(x / Lx) * Nx + y / Ly * Ny > 24, 0, 1))

grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(bottom))

#####
##### Setup atmospheric and oceanic forcing
#####

Uₐ = XFaceField(grid)
Vₐ = YFaceField(grid)

# Atmosphere - sea ice stress
τᵤ = Field(ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Uₐ)
τᵥ = Field(ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Vₐ)

# Initialize the wind stress (constant in time)
set!(Uₐ, (x, y) -> 𝓋ₐ)
compute!(τᵤ)
compute!(τᵥ)

#####
##### Ocean stress (a zero-velocity ocean with a drag)
#####
struct PrescribedOceanStress{FT}
    ρₒ :: FT
    Cᴰ :: FT
end

import ClimaSeaIce.SeaIceMomentumEquations: implicit_τx_coefficient, implicit_τy_coefficient

@inline function implicit_τx_coefficient(i, j, k, grid, τ::PrescribedOceanStress, clock, fields) 
    uᵢ = @inbounds fields.u[i, j, k]
    vᵢ = ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    
    return τ.ρₒ * τ.Cᴰ * sqrt(uᵢ^2 + vᵢ^2)
end

@inline function implicit_τy_coefficient(i, j, k, grid, τ::PrescribedOceanStress, clock, fields) 
    uᵢ = ℑxyᶠᶜᵃ(i, j, k, grid, fields.u)
    vᵢ = @inbounds fields.v[i, j, k]
    
    return τ.ρₒ * τ.Cᴰ * sqrt(uᵢ^2 + vᵢ^2)
end

τₒ = PrescribedOceanStress(1025.0, 5.5e-3)

#####
##### Numerical details
#####

# We use an elasto-visco-plastic rheology and WENO seventh order 
# for advection of h and ℵ
ice_dynamics = SeaIceMomentumEquation(grid; 
                                      coriolis = BetaPlane(latitude=60),
                                      solver = SplitExplicitSolver(substeps=120))
                                      
advection = WENO(; order = 7)

u_bcs = FieldBoundaryConditions(top = nothing, bottom = nothing,
                                north = ValueBoundaryCondition(0),
                                south = ValueBoundaryCondition(0))

#Define the model!
model = SeaIceModel(grid; 
                    top_external_stress = (u=τᵤ, v=τᵥ),
                    bottom_external_stress = (u=τₒ, v=τₒ), 
                    advection,
                    ice_dynamics = ice_dynamics,
                    boundary_conditions = (; u = u_bcs),
                    ice_thermodynamics = nothing)

# Initial height field with perturbations around 0.3 m
h₀(x, y) = 1.0

# We start with a concentration of ℵ = 1
set!(model, h = h₀)
set!(model, ℵ = 1)

#####
##### Setup the simulation
#####

# run the model for 2 days
simulation = Simulation(model, Δt = 2minutes, stop_time = 2days) 

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
    push!(htimeseries, deepcopy(Array(interior(h))))
    push!(ℵtimeseries, deepcopy(Array(interior(ℵ))))
    push!(utimeseries, deepcopy(Array(interior(u))))
    push!(vtimeseries, deepcopy(Array(interior(v))))
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

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e), max(trac): %.2f, %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., hmax, ℵmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))
simulation.callbacks[:save]     = Callback(accumulate_timeseries, IterationInterval(10))

run!(simulation)

# Visualize!
Nt = length(htimeseries)
iter = Observable(1)

hi = @lift(htimeseries[$iter][:, :, 1])
ℵi = @lift(ℵtimeseries[$iter][:, :, 1])
ui = @lift(utimeseries[$iter][:, :, 1])
vi = @lift(vtimeseries[$iter][:, :, 1])

fig = Figure(size = (1000, 500))
ax = Axis(fig[1, 1], title = "sea ice thickness")
heatmap!(ax, hi, colormap = :magma,         colorrange = (0.0, 2.0))

ax = Axis(fig[1, 2], title = "sea ice concentration")
heatmap!(ax, ℵi, colormap = Reverse(:deep), colorrange = (0.0, 1))

ax = Axis(fig[2, 1], title = "zonal velocity")
heatmap!(ax, ui, colorrange = (0, 0.12), colormap = :balance)

ax = Axis(fig[2, 2], title = "meridional velocity")
heatmap!(ax, vi, colorrange = (-0.025, 0.025), colormap = :bwr)

GLMakie.record(fig, "sea_ice_dynamics.mp4", 1:Nt, framerate = 50) do i
    iter[] = i
    @info "doing iter $i"
end
