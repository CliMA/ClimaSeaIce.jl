# # Sea ice advected by an atmospheric anticyclone
#
#
#
#
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using ClimaSeaIce
using Printf
using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies
using Oceananigans.Operators

# The experiment found in the paper: 
# Simulating Linear Kinematic Features in Viscous-Plastic Sea Ice Models 
# on Quadrilateral and Triangular Grids With Different Variable Staggering

arch = CPU()

L  = 512kilometers
𝓋ₒ = 0.01 # m / s maximum ocean speed
𝓋ₐ = 20.0 # m / s maximum atmospheric speed modifier
Cᴰ = 1.2e-3 # Atmosphere - sea ice drag coefficient
ρₐ = 1.3  # kg/m³

# 2 km domain
grid = RectilinearGrid(arch;
                       size = (256, 256), 
                          x = (0, L), 
                          y = (0, L), 
                       halo = (8, 8),
                   topology = (Bounded, Bounded, Flat))

#####
##### Ocean sea-ice stress
#####

using ClimaSeaIce.SeaIceMomentumEquations: SemiImplicitOceanSeaIceStress

# Constant ocean velocities corresponding to a cyclonic eddy
Uₒ = Field{Face, Face, Center}(grid)
Vₒ = Field{Face, Face, Center}(grid)

set!(Uₒ, (x, y) -> 𝓋ₒ * (2y - L) / L)
set!(Vₒ, (x, y) -> 𝓋ₒ * (L - 2x) / L)

Oceananigans.BoundaryConditions.fill_halo_regions!(Uₒ)
Oceananigans.BoundaryConditions.fill_halo_regions!(Vₒ)

τᵤₒ = τᵥₒ = SemiImplicitOceanSeaIceStress(Uₒ, Vₒ, 5.5e-3, 1026.0)

####
#### Atmosphere - sea ice stress 
####

Uₐ = Field{Face, Face, Center}(grid)
Vₐ = Field{Face, Face, Center}(grid)

# Atmosphere - sea ice stress
τᵤₐ = Field(ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Uₐ)
τᵥₐ = Field(ρₐ * Cᴰ * sqrt(Uₐ^2 + Vₐ^2) * Vₐ)

# Atmospheric velocities corresponding to an anticyclonic eddy moving north-east
@inline center(t) = 256kilometers + 51.2kilometers * t / 86400
@inline radius(x, y, t)  = sqrt((x - center(t))^2 + (y - center(t))^2)
@inline speed(x, y, t)   = exp(- radius(x, y, t) / 100kilometers) / 100 

@inline ua_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (cosd(72) * (x - center(t)) + sind(72) * (y - center(t))) / 1000
@inline va_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (cosd(72) * (y - center(t)) - sind(72) * (x - center(t))) / 1000

# Initialize the stress at time t = 0
set!(Uₐ, (x, y) -> ua_time(x, y, 0))
set!(Vₐ, (x, y) -> va_time(x, y, 0))
compute!(τᵤₐ)
compute!(τᵥₐ)

Oceananigans.BoundaryConditions.fill_halo_regions!(τᵤₐ)
Oceananigans.BoundaryConditions.fill_halo_regions!(τᵥₐ)

#####
##### Numerical details
#####

rheology = BrittleBinghamMaxwellRheology() 

# rheology =  ElastoViscoPlasticRheology(min_substeps=50, 
#                                        max_substeps=100,
#                                        minimum_plastic_stress=1e-10)

# We use an elasto-visco-plastic rheology and WENO seventh order 
# for advection of h and ℵ
momentum_equations = SeaIceMomentumEquation(grid; 
                                            top_momentum_stress = (u = τᵤₐ, v = τᵥₐ),
                                            bottom_momentum_stress = (u = τᵤₒ, v = τᵥₒ),
                                            coriolis = FPlane(f=1.56e-4),
                                            ocean_velocities = (u = Uₒ, v = Vₒ),
                                            rheology,
                                            solver = SplitExplicitSolver(substeps=150))

# Define the model!
model = SeaIceModel(grid; 
                    dynamics = momentum_equations,
                    ice_thermodynamics = nothing, # No thermodynamics here
                    advection = WENO(order=7),
                    timestepper = :QuasiAdamsBashforth2)

model.timestepper.χ = -0.5 # Euler forward

# Initial height field with perturbations around 0.3 m
h₀(x, y) = 0.3 + 0.005 * (sin(60 * x / 1000kilometers) + sin(30 * y / 1000kilometers))

# We start with a concentration of ℵ = 1
set!(model, h = h₀)
set!(model, ℵ = 1)

#####
##### Setup the simulation
#####

# run the model for 2 days
simulation = Simulation(model, Δt = 2minutes, stop_time = 2days) #, stop_iteration=1)

# Remember to evolve the wind stress field in time!
function compute_wind_stress(sim)
    time = sim.model.clock.time
    @inline ua(x, y) = ua_time(x, y, time)
    @inline va(x, y) = va_time(x, y, time)
    set!(Uₐ, ua)
    set!(Vₐ, va)

    compute!(τᵤₐ)
    compute!(τᵥₐ)

    Oceananigans.BoundaryConditions.fill_halo_regions!(τᵤₐ)
    Oceananigans.BoundaryConditions.fill_halo_regions!(τᵥₐ)
    
    return nothing
end

simulation.callbacks[:compute_stress] = Callback(compute_wind_stress, IterationInterval(1))

# Container to hold the data
htimeseries = []
ℵtimeseries = []
utimeseries = []
vtimeseries = []
σ₁₁timeseries = []
σ₁₂timeseries = []
σ₂₂timeseries = []
dtimeseries   = [] 

# Callback function to collect the data from the `sim`ulation
function accumulate_timeseries(sim)
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    u = sim.model.velocities.u
    v = sim.model.velocities.v

    σ₁₁ = sim.model.dynamics.auxiliary_fields.σ₁₁
    σ₁₂ = sim.model.dynamics.auxiliary_fields.σ₁₂
    σ₂₂ = sim.model.dynamics.auxiliary_fields.σ₂₂
    
    if haskey(sim.model.tracers, :d)
        d = sim.model.tracers.d
    else
        d = sim.model.ice_concentration
    end

    push!(htimeseries,   deepcopy(Array(interior(h))))
    push!(ℵtimeseries,   deepcopy(Array(interior(ℵ))))
    push!(utimeseries,   deepcopy(Array(interior(u))))
    push!(vtimeseries,   deepcopy(Array(interior(v))))
    push!(σ₁₁timeseries, deepcopy(Array(interior(σ₁₁))))
    push!(σ₁₂timeseries, deepcopy(Array(interior(σ₁₂))))
    push!(σ₂₂timeseries, deepcopy(Array(interior(σ₂₂))))
    push!(dtimeseries,   deepcopy(Array(interior(d  ))))
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

using GLMakie

# using JLD2
# jldsave("ice_anticyclone.jld2"; h=htimeseries, ℵ=ℵtimeseries, u=utimeseries, v=vtimeseries, σ₁₁=σ₁₁timeseries, σ₁₂=σ₁₂timeseries, σ₂₂=σ₂₂timeseries, d=dtimeseries)

# Visualize!
Nt = length(htimeseries)
iter = Observable(1)

hi   = @lift(htimeseries[$iter][:, :, 1])
ℵi   = @lift(ℵtimeseries[$iter][:, :, 1])
ui   = @lift(utimeseries[$iter][:, :, 1])
vi   = @lift(vtimeseries[$iter][:, :, 1])
σ₁₁i = @lift(σ₁₁timeseries[$iter][:, :, 1])
σ₁₂i = @lift(σ₁₂timeseries[$iter][:, :, 1])
σ₂₂i = @lift(σ₂₂timeseries[$iter][:, :, 1])
di   = @lift(dtimeseries[$iter][:, :, 1])

∂xu = @lift(∂x(utimeseries[$iter]))
∂yu = @lift(∂y(utimeseries[$iter]))
∂xv = @lift(∂x(vtimeseries[$iter]))
∂yv = @lift(∂y(vtimeseries[$iter]))

ϵ = @lift(interior(compute!(Field(sqrt(($∂xu + $∂yv)^2 + ($∂yu - $∂xv)^2))), :, :, 1))

fig = Figure()
ax = Axis(fig[1, 1], title = "sea ice thickness")
heatmap!(ax, hi, colormap = :magma,         colorrange = (0.23, 0.37))

ax = Axis(fig[1, 2], title = "sea ice concentration")
heatmap!(ax, ℵi, colormap = Reverse(:deep), colorrange = (0.75, 1))

ax = Axis(fig[2, 1], title = "damage")
heatmap!(ax, di, colorrange = (0.8, 1.0))

ax = Axis(fig[2, 2], title = "total deformation")
heatmap!(ax, ϵ, colorrange = (0, 1e-5), colormap = Reverse(:grays))

record(fig, "sea_ice_rheology.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
    @info "doing iter $i"
end

fig = Figure()
ax = Axis(fig[1, 1], title = "sigma 11")
heatmap!(ax, σ₁₁i)

ax = Axis(fig[1, 2], title = "sigma 22")
heatmap!(ax, σ₂₂i)

ax = Axis(fig[2, 1], title = "sigma 12")
heatmap!(ax, σ₁₂i)

ax = Axis(fig[2, 2], title = "damage")
heatmap!(ax, di, colorrange = (0.8, 1.0))

record(fig, "sea_ice_stress.mp4", 1:Nt-1, framerate = 8) do i
    iter[] = i
    @info "doing iter $i"
end
