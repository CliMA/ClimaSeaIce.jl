using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using Printf
using ClimaSeaIce.SlabSeaIceModels.SlabSeaIceDynamics.ExplicitRheologies: ElastoViscoPlasticRheology

# The experiment found in the paper: 
# Simulating Linear Kinematic Features in Viscous-Plastic Sea Ice Models 
# on Quadrilateral and Triangular Grids With Different Variable Staggering

L  = 512kilometers
𝓋ₒ = 0.01 # m / s maximum ocean speed
𝓋ₐ = 30.0 # m / s maximum atmospheric speed times `e`
Cᴰ = 1.2e-3 # Atmosphere - sea ice drag coefficient
ρₐ = 1.3 # kg/m³

# 2 km domain
grid = RectilinearGrid(size=(256, 256), 
                         x = (0, L), 
                         y = (0, L), 
                  topology = (Bounded, Bounded, Flat))

# Constant ocean velocities corresponding to a cyclonic eddy
Uₒ = XFaceField(grid)
Vₒ = YFaceField(grid)

set!(Uₒ, (x, y) -> 𝓋ₒ * (2y - L) / L)
set!(Vₒ, (x, y) -> 𝓋ₒ * (L - 2x) / L)

Uₐ = XFaceField(grid)
Vₐ = YFaceField(grid)

# Atmosphere - sea ice stress
τᵤ = Field(ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Uₐ)
τᵥ = Field(ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Vₐ)

# Atmospheric velocities corresponding to an anticyclonic eddy moving north-east
@inline center(t) = 256kilometers + 51.2kilometers * t / 86400
@inline radius(x, y, t)  = sqrt((x - center(t))^2 + (y - center(t))^2)
@inline speed(x, y, t)   = 1 / 100 * exp(- radius(x, y, t) / 100kilometers)

@inline ua_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (cosd(72) * (x - center(t)) + sind(72) * (y - center(t))) / 1000
@inline va_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (cosd(72) * (y - center(t)) - sind(72) * (x - center(t))) / 1000

# Initialize the stress at time t = 0
set!(Uₐ, (x, y) -> ua_time(x, y, 0))
set!(Vₐ, (x, y) -> va_time(x, y, 0))
compute!(τᵤ)
compute!(τᵥ)

# We use an elasto-visco-plastic rheology and WENO seventh order for advection
rheology  = ElastoViscoPlasticRheology(grid; substeps = 1000)
advection = Centered()

# Define the model!
model = SlabSeaIceModel(grid; 
                        top_u_stress = τᵤ,
                        top_v_stress = τᵥ,
                        ocean_velocities = (u = Uₒ, v = Vₒ),
                        rheology,
                        advection,
                        coriolis  = FPlane(f = 1e-4),
                        top_heat_boundary_condition=PrescribedTemperature(-10))

# Initial height field with perturbations around 0.3 m
h₀(x, y) = 0.3 + 0.005 * (sin(60 * x / 1000kilometers) + sin(30 * y / 1000kilometers))

# We start with a concentration of ℵ = 1
set!(model, h = h₀)
set!(model, ℵ = 1)

# run the model for 2 days
simulation = Simulation(model, Δt = 2minutes, stop_time = 10hours)

# Remember to evolve the wind stress field in time!
function compute_wind_stress(sim)
    time = sim.model.clock.time
    @inline ua(x, y) = ua_time(x, y, time)
    @inline va(x, y) = va_time(x, y, time)
    set!(Uₐ, ua)
    set!(Vₐ, va)

    compute!(τᵤ)
    compute!(τᵥ)

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

run!(simulation)