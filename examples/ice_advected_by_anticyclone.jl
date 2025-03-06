# # Sea ice advected by an atmospheric anticyclone
# 
#
#
#
using ClimaSeaIce
using ClimaSeaIce.SeaIceMomentumEquations
using ClimaSeaIce.Rheologies
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Printf
using CairoMakie

# The experiment found in the paper: 
# Simulating Linear Kinematic Features in Viscous-Plastic Sea Ice Models 
# on Quadrilateral and Triangular Grids With Different Variable Staggering

arch = CPU()

L  = 512kilometers
𝓋ₒ = 0.01 # m / s maximum ocean speed
𝓋ₐ = 30.0 # m / s maximum atmospheric speed modifier

# 2 km domain
grid = RectilinearGrid(arch;
                       size = (256, 256), 
                          x = (0, L), 
                          y = (0, L), 
                       halo = (7, 7),
                   topology = (Bounded, Bounded, Flat))

#####                   
##### Value boundary conditions for velocities
#####

u_bcs = FieldBoundaryConditions(north=ValueBoundaryCondition(0),
                                south=ValueBoundaryCondition(0))

v_bcs = FieldBoundaryConditions(west=ValueBoundaryCondition(0),
                                east=ValueBoundaryCondition(0))

#####
##### Ocean sea-ice stress
#####

# Constant ocean velocities corresponding to a cyclonic eddy
Uₒ = Field{Face, Face, Nothing}(grid)
Vₒ = Field{Face, Face, Nothing}(grid)

set!(Uₒ, (x, y) -> 𝓋ₒ * (2y - L) / L)
set!(Vₒ, (x, y) -> 𝓋ₒ * (L - 2x) / L)
fill_halo_regions!((Uₒ, Vₒ))

τₒ = SemiImplicitStress(uₑ=Uₒ, vₑ=Vₒ)

####
#### Atmosphere - sea ice stress 
####

Uₐ = Field{Face, Face, Nothing}(grid)
Vₐ = Field{Face, Face, Nothing}(grid)

τₐ = SemiImplicitStress(; uₑ=Uₐ, vₑ=Vₐ, ρₑ=1.3, Cᴰ=1.2e-3)

# Atmospheric velocities corresponding to an anticyclonic eddy moving north-east
@inline center(t) = 256kilometers + 51.2kilometers * t / 86400
@inline radius(x, y, t)  = sqrt((x - center(t))^2 + (y - center(t))^2)
@inline speed(x, y, t)   = 1 / 100 * exp(- radius(x, y, t) / 100kilometers)

@inline ua_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (cosd(72) * (x - center(t)) + sind(72) * (y - center(t))) / 1000
@inline va_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (cosd(72) * (y - center(t)) - sind(72) * (x - center(t))) / 1000

# Initialize the stress at time t = 0
set!(Uₐ, (x, y) -> ua_time(x, y, 0))
set!(Vₐ, (x, y) -> va_time(x, y, 0))

fill_halo_regions!((Uₐ, Vₐ))

#####
##### Numerical details
#####

# We use an elasto-visco-plastic rheology and WENO seventh order 
# for advection of h and ℵ

momentum_equations = SeaIceMomentumEquation(grid; 
                                            top_momentum_stress = τₐ,
                                            bottom_momentum_stress = τₒ,
                                            coriolis = FPlane(f=1e-4),
                                            ocean_velocities = (u = Uₒ, v = Vₒ),
                                            rheology = BrittleBinghamMaxwellRheology(),
                                            solver   = SplitExplicitSolver(substeps=150))

# Define the model!

model = SeaIceModel(grid; 
                    dynamics = momentum_equations,
                    ice_thermodynamics = nothing, # No thermodynamics here
                    advection = WENO(order=7),
                    tracers = :d,
                    boundary_conditions = (u=u_bcs, v=v_bcs))

# We start with a concentration of ℵ = 1 and an 
# initial height field with perturbations around 0.3 m

h₀(x, y) = 0.3 + 0.005 * (sin(60 * x / 1000kilometers) + sin(30 * y / 1000kilometers))

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

    fill_halo_regions!((Uₐ, Vₐ))
    
    return nothing
end

simulation.callbacks[:top_stress] = Callback(compute_wind_stress, IterationInterval(1))

h = model.ice_thickness
ℵ = model.ice_concentration
u, v = model.velocities

outputs = (; h, u, v, ℵ)

simulation.output_writers[:sea_ice] = JLD2OutputWriter(model, outputs;
                                                       filename = "sea_ice_advected_by_anticyclone.jld2", 
                                                       schedule = IterationInterval(5),
                                                       overwrite_existing = true)

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

run!(simulation)

using CairoMakie

htimeseries = FieldTimeSeries("sea_ice_advected_by_anticyclone.jld2", "h")
utimeseries = FieldTimeSeries("sea_ice_advected_by_anticyclone.jld2", "u")
vtimeseries = FieldTimeSeries("sea_ice_advected_by_anticyclone.jld2", "v")
ℵtimeseries = FieldTimeSeries("sea_ice_advected_by_anticyclone.jld2", "ℵ")

# Visualize!
Nt = length(htimeseries)
iter = Observable(1)

hi = @lift(htimeseries[$iter])
ℵi = @lift(ℵtimeseries[$iter])
ui = @lift(utimeseries[$iter])
vi = @lift(vtimeseries[$iter])

fig = Figure()
ax = Axis(fig[1, 1], title = "sea ice thickness")
heatmap!(ax, hi, colormap = :magma, colorrange = (0.23, 0.37))

ax = Axis(fig[1, 2], title = "sea ice concentration")
heatmap!(ax, ℵi, colormap = Reverse(:deep), colorrange = (0.9, 1))

ax = Axis(fig[2, 1], title = "zonal velocity")
heatmap!(ax, ui, colorrange = (-0.1, 0.1))

ax = Axis(fig[2, 2], title = "meridional velocity")
heatmap!(ax, vi, colorrange = (-0.1, 0.1))

CairoMakie.record(fig, "sea_ice_advected_by_anticyclone.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
    @info "doing iter $i"
end
nothing #hide

# ![](sea_ice_advected_by_anticyclone.mp4)
