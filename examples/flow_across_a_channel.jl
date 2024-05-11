using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using Printf
using ClimaSeaIce.SlabSeaIceModels.SlabSeaIceDynamics.ExplicitRheologies: ElastoViscoPlasticRheology

Lx = 200kilometers
Ly = 100kilometers

grid = RectilinearGrid(size=(200, 200), x = (0, Lx), y = (0, Ly), topology = (Periodic, Bounded, Flat))

bottom(x, y) = (- 10000 < x - Lx < 10000) & (y < Ly / 3) 

# Have an immersed boundary in the middle of the domain
grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(bottom))

# the atmosphere has a constant velocity pushing on the ice of 5 m/s
τₐ = XFaceField(grid)
uₐ = 5 # m/s
Cᴰ = 1e-3 
ρₐ = 1.225 # kg/m³
set!(τₐ, ρₐ * 1.225 * uₐ^2)

Oceananigans.fill_halo_regions!(τₐ)

rheology = ElastoViscoPlasticRheology(grid; substeps = 200)

model = SlabSeaIceModel(grid; 
                        top_u_stress = τₐ,
                        rheology,
                        advection = WENO(; order = 7),
                        top_heat_boundary_condition=PrescribedTemperature(-10))

set!(model, h = 1)
set!(model, ℵ = 1)

# run the model for 3 days
simulation = Simulation(model, Δt = 1minutes)

## Container to hold the data
htimeseries = []
ℵtimeseries = []
utimeseries = []
vtimeseries = []

## Callback function to collect the data from the `sim`ulation
function accumulate_timeseries(sim)
    h = sim.model.ice_thickness
    ℵ = sim.model.concentration
    u = sim.model.velocities.u
    v = sim.model.velocities.v
    push!(htimeseries, interior(h))
    push!(ℵtimeseries, interior(ℵ))
    push!(utimeseries, interior(u))
    push!(vtimeseries, interior(v))
end

wall_time = [time_ns()]

function progress(sim) 
    h = sim.model.ice_thickness
    ℵ = sim.model.concentration
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
simulation.callbacks[:save]     = Callback(accumulate_timeseries, IterationInterval(5))
