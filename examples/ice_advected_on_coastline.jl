using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics
using ClimaSeaIce.Rheologies
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Printf
using CairoMakie

# A solid block of ice moving against a triangular coastline in a periodic channel

Lx = 512kilometers
Ly = 256kilometers
Nx = 256
Ny = 256

y_max = Ly / 2

arch = CPU()

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

# Atmosphere - sea ice stress
Ua = XFaceField(grid)
τᵤ = Field(- ρₐ * Cᴰ * Ua^2)
set!(Ua, 𝓋ₐ)
compute!(τᵤ)
Oceananigans.BoundaryConditions.fill_halo_regions!(τᵤ)
τᵥ = 0.0

#####
##### Ocean stress (a zero-velocity ocean with a drag)
#####

τₒ = SemiImplicitStress()

#####
##### Numerical details
#####

# We use an elasto-visco-plastic rheology and WENO seventh order 
# for advection of h and ℵ
dynamics = SeaIceMomentumEquation(grid; 
                                  top_momentum_stress = (u=τᵤ, v=τᵥ),
                                  bottom_momentum_stress = τₒ, 
                                  rheology = ElastoViscoPlasticRheology(),
                                  solver = SplitExplicitSolver(substeps=150))

@inline immersed_u_drag(i, j, k, grid, clock, fields, D) = @inbounds - D * fields.u[i, j, k]
@inline immersed_v_drag(i, j, k, grid, clock, fields, D) = @inbounds - D * fields.v[i, j, k]

immersed_u_bc = FluxBoundaryCondition(immersed_u_drag, discrete_form=true, parameters=3e-2)
immersed_v_bc = FluxBoundaryCondition(immersed_v_drag, discrete_form=true, parameters=3e-2)

immersed_u_bc = ImmersedBoundaryCondition(top=nothing, bottom=nothing, west=nothing, east=nothing, 
                                          south=immersed_u_bc, north=immersed_u_bc)

immersed_v_bc = ImmersedBoundaryCondition(top=nothing, bottom=nothing, south=nothing, north=nothing,
                                          west=immersed_v_bc, east=immersed_v_bc)

u_bcs = FieldBoundaryConditions(grid, (Face, Center, Nothing); 
                                north = ValueBoundaryCondition(0),
                                south = ValueBoundaryCondition(0),
                                immersed = immersed_u_bc)

v_bcs = FieldBoundaryConditions(grid, (Center, Face, Nothing); 
                                immersed = immersed_v_bc)

#Define the model! 
model = SeaIceModel(grid; 
                    advection = WENO(order=7),
                    dynamics = dynamics,
                    boundary_conditions = (; u=u_bcs, v=v_bcs),
                    ice_thermodynamics = nothing)

# We start with a concentration of ℵ = 1 everywhere
set!(model, h = 1)
set!(model, ℵ = 1)

#####
##### Setup the simulation
#####

# run the model for 10 day
simulation = Simulation(model, Δt = 2minutes, stop_time=2days) 

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

wall_time = Ref(time_ns())

function progress(sim) 
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    u = sim.model.velocities.u
    v = sim.model.velocities.v

    hmax = maximum(interior(h))
    ℵmin = minimum(interior(ℵ))
    umax = maximum(interior(u)), maximum(interior(v))
    step_time = 1e-9 * (time_ns() - wall_time[])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e), max(trac): %.2f, %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., hmax, ℵmin, prettytime(step_time))

     wall_time[] = time_ns()
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

CairoMakie.record(fig, "sea_ice_advected_on_coastline.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
    @info "doing iter $i"
end
nothing #hide

# ![](sea_ice_advected_on_coastline.mp4)
