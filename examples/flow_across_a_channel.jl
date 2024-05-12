using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using Printf
using ClimaSeaIce.SlabSeaIceModels.SlabSeaIceDynamics.ExplicitRheologies: ElastoViscoPlasticRheology

Lx = 200kilometers
Ly = 200kilometers
σ² = 0.2

grid = RectilinearGrid(size=(100, 100), x = (0, Lx), y = (0, Ly), topology = (Bounded, Bounded, Flat))

# the atmosphere has a constant velocity pushing on the ice
# The ocean will have zero velocity and therefore "slow down" the ice
τᵤ = XFaceField(grid)
τᵥ = YFaceField(grid)
Uₐ = XFaceField(grid)
Vₐ = YFaceField(grid)

norm_x(x) = 2 * (x - Lx / 2) / Lx # from -1 to 1
norm_y(y) = 2 * (y - Ly / 2) / Ly # from -1 to 1

@inline η(x, y) = exp(-(norm_x(x)^2 + norm_y(y)^2) ./ σ²) 

""" background and initial zonal velocity """
@inline function U̅(x, y)
    f = 1e-4
    g = 10
    return - 5 * g * 2 * norm_y(y) * η(x, y) / f / Ly * 2
end

""" background and initial meridional velocity """
@inline function V̅(x, y)
    f  = 1e-4
    g  = 10
    return 5 * g * 2 * norm_x(x) * η(x, y) / f / Lx * 2
end

set!(Uₐ, U̅)
set!(Vₐ, V̅)

Cᴰ = 1e-3 
ρₐ = 1.225 # kg/m³
set!(τᵤ, ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Uₐ)
set!(τᵥ, ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Vₐ)

Oceananigans.fill_halo_regions!((τᵤ, τᵥ))

rheology = ElastoViscoPlasticRheology(grid; substeps = 300)

model = SlabSeaIceModel(grid; 
                        top_u_stress = τᵤ,
                        top_v_stress = τᵥ,
                        rheology,
                        advection = WENO(; order = 7),
                        coriolis  = FPlane(f = 1e-4),
                        top_heat_boundary_condition=PrescribedTemperature(-10))

set!(model, h = 1)
set!(model, ℵ = 1)

# run the model for 3 days
simulation = Simulation(model, Δt = 1minutes)

# Container to hold the data
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

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))
simulation.callbacks[:save]     = Callback(accumulate_timeseries, IterationInterval(5))
