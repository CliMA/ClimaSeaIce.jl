using Oceananigans

# This files sets up a model that resembles the Antarctic Slope Current (ASC) model in the
# 2022 paper by Ian Eisenman

arch = CPU()

g_Earth = 9.80665

Lx = 400000
Ly = 450000
Lz = 4000

Nx = 400
Ny = 450
Nz = 70 # TODO: modify spacing if needed, 10 m at surface, 100m at seafloor

#
# Setting up the grid and bathymetry:
#
# Using a linear slope that approximates the min and max spacings in the paper:
init_adjustment = (0.5*(100-10)/(Nz-1)) + (10 - (100-10)/(Nz-1))
linear_slope(k) = (0.5*(100-10)/(Nz-1))*(k)^2 + (10 - (100-10)/(Nz-1))*(k) - init_adjustment

spacing_adjustment      = Lz / linear_slope(Nz+1)
linear_slope_z_faces(k) = -spacing_adjustment * linear_slope(k)
# Can reverse this to get grid from -4000 to 0 later

underlying_grid = RectilinearGrid(arch,
                                  size = (Nx, Ny, Nz),
                                  topology = (Periodic, Bounded, Bounded),
                                  x = (-Lx/2, Lx/2),
                                  y = (0, Ly),
                                  z = linear_slope_z_faces)

println(underlying_grid)

# We want the underwater slope to provide a depth of 500 m at y = 0 and the full 4 km by y =200. It follows
# a hyperbolic tangent curve with its center at y ~= 150 at a depth of ~ -4000 + (4000 - 500)/2 = -2250 m
underwater_slope(x, y) = -1750tanh(0.00004*(y - 150000)) - 2250
# TODO: add 50km troughs to the slope, dependent on x and y. Use appendix C to add detail here

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(underwater_slope))

println(grid)

# TODO: add underwater slope at y = 0 with troughs

#
# Running the Hydrostatic model:
#

# Forcings:
u_forcing(x, y, z, t) = exp(z) * cos(x) * sin(t)

# Boundary Conditions:
no_slip_bc   = ValueBoundaryCondition(0.0)
free_slip_bc = FluxBoundaryCondition(nothing)

free_slip_surface_bcs = FieldBoundaryConditions(no_slip_bc, top=FluxBoundaryCondition(nothing))
no_slip_field_bcs     = FieldBoundaryConditions(no_slip_bc)


# Assuming no particles or biogeochemistry
model = HydrostaticFreeSurfaceModel(; grid,
                                         clock = Clock{eltype(grid)}(0, 0, 1),
                            momentum_advection = CenteredSecondOrder(),
                              tracer_advection = CenteredSecondOrder(),
                                      buoyancy = SeawaterBuoyancy(eltype(grid)),
                                      coriolis = nothing,
                                  free_surface = ImplicitFreeSurface(gravitational_acceleration=g_Earth),
                                       #forcing = (,), # NamedTuple(),
                                       closure = nothing,
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
