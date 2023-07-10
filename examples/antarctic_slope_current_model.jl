
#using CairoMakie
using Oceananigans
using Oceananigans.Units: minute, minutes, hour
using SeawaterPolynomials.TEOS10

# This file sets up a model that resembles the Antarctic Slope Current (ASC) model in the
# 2022 paper by Ian Eisenman

arch = GPU()

g_Earth = 9.80665

Lx = 400kilometers
Ly = 450kilometers
Lz = 4000

Nx = 400
Ny = 450
Nz = 70 # TODO: modify spacing if needed, 10 m at surface, 100m at seafloor

sponge_width = 20000
#=
#
# Setting up the grid and bathymetry:
#
# Using a linear slope that approximates the min and max spacings in the paper:
init_adjustment = (0.5*(100-10)/(Nz-1)) + (10 - (100-10)/(Nz-1))
linear_slope(k) = (0.5*(100-10)/(Nz-1))*(k)^2 + (10 - (100-10)/(Nz-1))*(k) - init_adjustment

spacing_adjustment      = Lz / linear_slope(Nz+1)
linear_slope_z_faces(k) = -spacing_adjustment * linear_slope(k)
# Can reverse this to get grid from -4000 to 0 later
=#
refinement = 20.0 # controls spacing near surface (higher means finer spaced), 12.0
stretching = 3  # controls rate of stretching at bottom, 3

# Normalized height ranging from 0 to 1
h(k) = (k - 1) / Nz

# Linear near-surface generator
ζ₀(k) = 1 + (h(k) - 1) / refinement

# Bottom-intensified stretching function
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

# Generating function
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

underlying_grid = RectilinearGrid(arch,
                                  size = (Nx, Ny, Nz),
                                  topology = (Periodic, Bounded, Bounded),
                                  x = (-Lx/2, Lx/2),
                                  y = (0, Ly),
                                  z = z_faces)

@show underlying_grid

# We want the underwater slope to provide a depth of 500 m at y = 0 and the full 4 km by y =200. It follows
# a hyperbolic tangent curve with its center at y ~= 150 at a depth of ~ -4000 + (4000 - 500)/2 = -2250 m
y₀ = 150kilometers
Δy = 25kilometers

slope_depth = 500
basin_depth = 4000

""" Varies smoothly between 0 when y << y₀ to 1 when y >> y₀ + Δy, over a region of width Δy. """
step(y, y₀, Δy) = (1 + tanh((y - y₀) / Δy)) / 2

underwater_slope(x, y) = -slope_depth + (slope_depth - basin_depth) * tanh((y - y₀) / Δy) 
# TODO: add 50km troughs to the slope, dependent on x and y. Use appendix C to add detail here

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(underwater_slope))

@show grid

# TODO: add underwater slope at y = 0 with troughs

#
# Creating the Hydrostatic model:
#

# Forcings:
u_forcing(x, y, z, t) = exp(z) * cos(x) * sin(t) # Not actual forcing, just example

# TODO: need to use forcings to enact the sponge layers. Their support is within 20000 m
# of the N/S boundaries, they impose a cross-slope buoyancy gradient, and the relaxation
# tiem scales decrease linearly with distance from the interior termination of the sponge layers

# We'll use Relaxation() to impose a sponge layer forcing on velocity, temperature, and salinity
damping_rate       = 1 / 43200 # Relaxation time scale in seconds, need to make this decrease linearly toward outermost boundary
south_mask         = GaussianMask{:y}(center=0, width=sponge_width)
north_mask         = GaussianMask{:y}(center=Ly, width=sponge_width)
south_sponge_layer = Relaxation(; rate=damping_rate, mask=south_mask)
north_sponge_layer = Relaxation(; rate=damping_rate, mask=north_mask)
sponge_layers      = (south_sponge_layer, north_sponge_layer)
# TODO: compose north_mask and south_mask together into one sponge layer, OR compose north/south sponge layers


# Boundary Conditions:
no_slip_bc   = ValueBoundaryCondition(0.0)
free_slip_bc = FluxBoundaryCondition(nothing)

free_slip_surface_bcs = FieldBoundaryConditions(no_slip_bc, top=FluxBoundaryCondition(nothing))
no_slip_field_bcs     = FieldBoundaryConditions(no_slip_bc)

# Buoyancy Equations of State - we want high order polynomials, so we'll use TEOS-10
eos = TEOS10EquationOfState()

# Coriolis Effect, using basic f-plane with precribed reference Coriolis parameter
coriolis = FPlane(f=-1.3e14)

# Diffusivities as part of closure
# TODO: make sure this works for biharmonic diffusivities as the horizontal, 
horizontal_closure = HorizontalScalarDiffusivity(ν=0.1, κ=0.1)
vertical_closure   = VerticalScalarDiffusivity(ν=3e-4, κ=1e-5)

# Assuming no particles or biogeochemistry
model = HydrostaticFreeSurfaceModel(; grid,
                                         clock = Clock{eltype(grid)}(0, 0, 1),
                            momentum_advection = WENO(),
                              tracer_advection = WENO(),
                                      buoyancy = SeawaterBuoyancy(equation_of_state=eos, gravitational_acceleration=g_Earth),
                                      coriolis = coriolis,
                                  free_surface = ImplicitFreeSurface(gravitational_acceleration=g_Earth),
                                       forcing = (u=sponge_layers, v=sponge_layers, w=sponge_layers, T=sponge_layers, S=sponge_layers), # NamedTuple()
                                       closure = CATKEVerticalDiffusivity(),
                           boundary_conditions = (u=free_slip_surface_bcs, v=free_slip_surface_bcs, w=no_slip_field_bcs), # NamedTuple(),
                                       tracers = (:T, :S),
                                    velocities = nothing,
                                      pressure = nothing,
                            diffusivity_fields = nothing,
                              #auxiliary_fields = nothing, # NamedTuple(),
)

println(model)

println("u boundary conditions:")
println(model.velocities.u.boundary_conditions)
println("v boundary conditions:")
println(model.velocities.v.boundary_conditions)
println("w boundary conditions:")
println(model.velocities.w.boundary_conditions)

#
# Now create a simulation and run the model
#
simulation = Simulation(model; Δt=100.0, stop_time=600minutes)

# Create a NamedTuple

filename = "asc_model_run"

simulation.output_writers[:slices] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                     filename = filename * ".jld2",
                     indices = (:, grid.Ny/2, :),
                     schedule = TimeInterval(1minute),
                     overwrite_existing = true)

run!(simulation)

@show simulation
