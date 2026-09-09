# # Ice advected past an obstacle
#
# This example simulates a solid sheet of ice driven by wind past a rhomboidal island in a doubly
# periodic domain. The island is represented by an immersed boundary, and the ice sticks to it: the
# island's faces impose no slip, so the ice piles up on the upwind side and forms a wake downwind.
# This example demonstrates how to:
#
#   * set up a two-dimensional sea ice model with immersed boundaries,
#   * impose no slip on an immersed boundary with a zero-valued boundary condition,
#   * prescribe atmospheric and oceanic stresses,
#   * use elasto-visco-plastic rheology with split-explicit time stepping,
#   * visualize the evolution of ice thickness, concentration, and velocity.
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
# We set up a doubly periodic two-dimensional domain holding a single island:

using ClimaSeaIce
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Printf
using CairoMakie

Lx = 512kilometers
Ly = 256kilometers
Nx = 256
Ny = 128

arch = CPU()

# ## Grid configuration
#
# We create a rectilinear grid that is periodic in both ``x`` and ``y``, so the only boundary the
# ice can feel is the island itself:

grid = RectilinearGrid(arch; size = (Nx, Ny),
                                x = (-Lx/2, Lx/2),
                                y = (0, Ly),
                             halo = (4, 4),
                         topology = (Periodic, Periodic, Flat))

# The island is a rhomboid centred in the domain, described by the set of points whose distance from
# the centre, measured in the ``L¹`` norm, is less than one:

island_half_width = 48kilometers
island_half_height = 32kilometers
island_centre = (0, Ly/2)

@inline island(x, y) = abs(x - island_centre[1]) / island_half_width +
                       abs(y - island_centre[2]) / island_half_height < 1

grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(island))

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

# The island holds the ice at rest along its faces. As on a domain boundary, no slip is requested
# with a zero-valued boundary condition — here in the `immersed` slot. The rheology reads it and
# doubles the wall contribution to the shear strain rate, which is the gradient the reflected
# velocity inside the island would produce. Leaving `immersed` unset would instead let the ice slide
# freely along the island.

u_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing);
                                immersed = ValueBoundaryCondition(0))

v_bcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
                                immersed = ValueBoundaryCondition(0))

# We define the model with WENO advection and no thermodynamics:

model = SeaIceModel(grid;
                    advection = WENO(order=7),
                    dynamics,
                    boundary_conditions = (; u=u_bcs, v=v_bcs),
                    ice_thermodynamics = nothing)

# We initialize the model with uniform ice thickness and concentration:

set!(model, h = 1)
set!(model, ℵ = 1)

# ## Running the simulation
#
# We run the model for 3 days with a 5-minute time step:

simulation = Simulation(model, Δt=5minutes, stop_time=3days)

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

simulation.callbacks[:save] = Callback(accumulate_timeseries, IterationInterval(10))

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

CairoMakie.record(fig, "sea_ice_advected_past_an_obstacle.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
    @info "Rendering frame $i"
end
nothing #hide

# ![](sea_ice_advected_past_an_obstacle.mp4)
#
# The animation shows the ice sheet deforming as it is driven past the island. Ice converges and
# thickens on the upwind face, and the no-slip condition drags it to a halt along the island's flanks,
# leaving a slow wake that persists downwind.
