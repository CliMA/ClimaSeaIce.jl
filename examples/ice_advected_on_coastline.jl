# # Ice advected on coastline example
#
# This example simulates a solid block of ice moving against a triangular coastline
# in a periodic channel. The ice is driven by atmospheric winds and interacts with
# the coastline through immersed boundaries. This example demonstrates
#
#   * How to set up a two-dimensional sea ice model with immersed boundaries.
#   * How to prescribe atmospheric and oceanic stresses.
#   * How to use elasto-visco-plastic rheology with split-explicit time stepping.
#   * How to visualize the evolution of ice thickness, concentration, and velocity.
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
# We set up a two-dimensional periodic channel with a triangular coastline:

using ClimaSeaIce
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Printf
using CairoMakie

Lx = 512kilometers
Ly = 256kilometers
Nx = 256
Ny = 256

y_max = Ly / 2

arch = CPU()

# ## Grid configuration
#
# We create a rectilinear grid with periodic boundaries in x and bounded boundaries in y:

grid = RectilinearGrid(arch; size = (Nx, Ny), 
                                x = (-Lx/2, Lx/2), 
                                y = (0, Ly), 
                             halo = (4, 4),
                         topology = (Periodic, Bounded, Flat))

# We define a triangular coastline using an immersed boundary:

bottom(x, y) = ifelse(y > y_max, 0, 
               ifelse(abs(x / Lx) * Nx + y / Ly * Ny > 24, 0, 1))

grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(bottom))

# ## Atmospheric and oceanic forcing
#
# We set up atmospheric wind stress. The atmosphere has a constant wind speed:

𝓋ₐ = 10.0   # m s⁻¹
Cᴰ = 1.2e-3 # atmosphere-sea ice drag coefficient
ρₐ = 1.3    # kg m⁻³

# We create a field for the atmospheric wind and compute the stress:

Ua = XFaceField(grid)
τᵤ = Field(- ρₐ * Cᴰ * Ua^2)
set!(Ua, 𝓋ₐ)
compute!(τᵤ)
Oceananigans.BoundaryConditions.fill_halo_regions!(τᵤ)
τᵥ = 0.0

# The ocean stress is represented by a semi-implicit stress with zero ocean velocity:

τₒ = SemiImplicitStress()

# ## Model configuration
#
# We use an elasto-visco-plastic rheology and WENO seventh order for advection
# of ice thickness and concentration:

dynamics = SeaIceMomentumEquation(grid; 
                                  top_momentum_stress = (u=τᵤ, v=τᵥ),
                                  bottom_momentum_stress = τₒ, 
                                  rheology = ElastoViscoPlasticRheology(),
                                  solver = SplitExplicitSolver(substeps=150))

# We set boundary conditions for the velocity field:

u_bcs = FieldBoundaryConditions(top = nothing, bottom = nothing,
                                north = ValueBoundaryCondition(0),
                                south = ValueBoundaryCondition(0))

# We define the model with WENO advection and no thermodynamics:

model = SeaIceModel(grid; 
                    advection = WENO(order=7),
                    dynamics = dynamics,
                    boundary_conditions = (; u=u_bcs),
                    ice_thermodynamics = nothing)

# We initialize the model with uniform ice thickness and concentration:

set!(model, h = 1)
set!(model, ℵ = 1)

# ## Running the simulation
#
# We run the model for 2 days with a 2-minute time step:

simulation = Simulation(model, Δt = 2minutes, stop_time=2days) 

# ## Collecting data
#
# We set up containers to hold the time series data:

htimeseries = []
ℵtimeseries = []
utimeseries = []
vtimeseries = []

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

simulation.callbacks[:save]     = Callback(accumulate_timeseries, IterationInterval(10))

run!(simulation)

# ## Visualizing the results
#
# We create an animation of the ice thickness, concentration, and velocity fields:

Nt = length(htimeseries)
iter = Observable(1)

hi = @lift(htimeseries[$iter][:, :, 1])
ℵi = @lift(ℵtimeseries[$iter][:, :, 1])
ui = @lift(utimeseries[$iter][:, :, 1])
vi = @lift(vtimeseries[$iter][:, :, 1])

fig = Figure(size = (1000, 500))
ax = Axis(fig[1, 1], title = "Sea ice thickness (m)")
heatmap!(ax, hi, colormap = :magma, colorrange = (0.0, 2.0))

ax = Axis(fig[1, 2], title = "Sea ice concentration")
heatmap!(ax, ℵi, colormap = Reverse(:deep), colorrange = (0.0, 1))

ax = Axis(fig[2, 1], title = "Zonal velocity (m s⁻¹)")
heatmap!(ax, ui, colorrange = (0, 0.12), colormap = :balance)

ax = Axis(fig[2, 2], title = "Meridional velocity (m s⁻¹)")
heatmap!(ax, vi, colorrange = (-0.025, 0.025), colormap = :bwr)

CairoMakie.record(fig, "sea_ice_advected_on_coastline.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
    @info "Rendering frame $i"
end
nothing #hide

# ![](sea_ice_advected_on_coastline.mp4)
#
# The animation shows how the ice block moves and deforms as it interacts with
# the triangular coastline. The ice accumulates and thickens near the coast due
# to the no-slip boundary conditions.
