
#using CairoMakie
using Oceananigans
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity # FIND submodule with CATKE
using Oceananigans.Units: minute, minutes, hour, kilometers
using SeawaterPolynomials.TEOS10

using GLMakie
using Printf

# This file sets up a model that resembles the Antarctic Slope Current (ASC) model in the
# 2022 paper by Si, Stewart, and Eisenman

arch = CPU()

g_Earth = 9.80665

Lx = 400kilometers
Ly = 450kilometers
Lz = 4000

Nx = 64 #400
Ny = 64 #450
Nz = 32 #70 # TODO: modify spacing if needed, 10 m at surface, 100m at seafloor

sponge_width = 20kilometers
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
                                         z = z_faces,
                                      halo = (4, 4, 4))

@show underlying_grid

# We want the underwater slope to provide a depth of 500 m at y = 0 and the full 4 km by y =200. It follows
# a hyperbolic tangent curve with its center at y ~= 150 at a depth of ~ -4000 + (4000 - 500)/2 = -2250 m
y₀ = 150kilometers
Δy = 25kilometers

slope_depth = 500
basin_depth = 4000

""" Varies smoothly between 0 when y << y₀ to 1 when y >> y₀ + Δy, over a region of width Δy. """
step(y, y₀, Δy) = (1 + tanh((y - y₀) / Δy)) / 2

#underwater_slope(x, y) = -(slope_depth) + (slope_depth - basin_depth) * tanh((y - y₀) / Δy)
underwater_slope(x, y) = (-basin_depth - slope_depth + (slope_depth - basin_depth) * tanh((y - y₀) / Δy)) / 2
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

#= Currently not used
# Boundary Conditions:
# Eventually want open boundary conditions with tides, but can add that in later
no_slip_bc   = ValueBoundaryCondition(0.0)
free_slip_bc = FluxBoundaryCondition(nothing)

free_slip_surface_bcs = FieldBoundaryConditions(no_slip_bc, top=FluxBoundaryCondition(nothing))
no_slip_field_bcs     = FieldBoundaryConditions(no_slip_bc)
=#
# Buoyancy Equations of State - we want high order polynomials, so we'll use TEOS-10
eos = TEOS10EquationOfState() # can compare to linear EOS later (linear not recommended for polar regions)

# Coriolis Effect, using basic f-plane with precribed reference Coriolis parameter
coriolis = FPlane(latitude=-60)

# Diffusivities as part of closure
# TODO: make sure this works for biharmonic diffusivities as the horizontal, 
horizontal_closure = HorizontalScalarDiffusivity(ν=0.1, κ=0.1)
vertical_closure   = VerticalScalarDiffusivity(ν=3e-4, κ=1e-5)

# Assuming no particles or biogeochemistry
model = HydrostaticFreeSurfaceModel(;     grid = grid,
                            momentum_advection = WENO(),
                              tracer_advection = WENO(),
                                      buoyancy = SeawaterBuoyancy(equation_of_state=eos, gravitational_acceleration=g_Earth),
                                      coriolis = coriolis,
                                  free_surface = ImplicitFreeSurface(gravitational_acceleration=g_Earth),
                                       #forcing = (u=sponge_layers, v=sponge_layers, w=sponge_layers, T=sponge_layers, S=sponge_layers), # NamedTuple()
                                       closure = CATKEVerticalDiffusivity(),
                           #boundary_conditions = (u=free_slip_surface_bcs, v=free_slip_surface_bcs), # NamedTuple(),
                                       tracers = (:T, :S, :e)
)

# INITIAL CONDITION FOR TEMPERATURE: a baroclinically unstable temperature distribution
"""
    ramp(y, Δy)

Linear ramp from 0 to 1 between Ly/2 - Δy/2 and Ly/2 + Δy/2.

For example:
```
                y < Ly/2 - Δy/2 => ramp = 0
  Ly/2 - Δy/2 < y < Ly/2 + Δy/2 => ramp = y / Δy
                y > Ly/2 + Δy/2 => ramp = 1
```
"""
ramp(y, Δy) = min(max(0, (y - Ly/2)/Δy + 1/2), 1)

N² = 1e-5 # [s⁻²] vertical stratification
M² = 1e-7 # [s⁻²] meridional temperature gradient

Δy = 100kilometers # width of the region of the front
ΔT = Δy * M²       # temperature jump associated with the front
ϵT = 1e-2 * ΔT     # noise amplitude

Tᵢ(x, y, z) = N² * z + ΔT * ramp(y, Δy) + ϵT * randn()

set!(model, T=Tᵢ)

#=
using GLMakie

# Build coordinates with units of kilometers
x, y, z = 1e-3 .* nodes(grid, (Center(), Center(), Center()))

T = model.tracers.T

fig, ax, hm = heatmap(y, z, interior(T)[1, :, :],
                      colormap=:deep,
                      axis = (xlabel = "y [km]",
                              ylabel = "z [km]",
                              title = "T(x=0, y, z, t=0)",
                              titlesize = 24))

Colorbar(fig[1, 2], hm, label = "[m s⁻²]")

current_figure() # hide
fig
=#

# TODO: impose a temperature restoring at the surface that corresponds to Si/Stewart to help get baroclinic eddies
# Add initial condition to get turbulence running (should be doable on a laptop), then:
# Buoyancy gradient imposed at North/South boundary, surface temp as freezing temperature imposed

@show model

#
# Now create a simulation and run the model
#
simulation = Simulation(model; Δt=100.0, stop_iteration=100)

# Create a NamedTuple

filename = "asc_model_run"

simulation.output_writers[:slices] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                     filename = filename * ".jld2",
                     indices = (:, :, grid.Nz),
                     schedule = IterationInterval(1),
                     overwrite_existing = true)

run!(simulation)

@show simulation

#
# Make a figure and plot it
#

filepath = filename * ".jld2"

time_series = (u = FieldTimeSeries(filepath, "u"),
               v = FieldTimeSeries(filepath, "v"),
               w = FieldTimeSeries(filepath, "w"),
               T = FieldTimeSeries(filepath, "T"),
               S = FieldTimeSeries(filepath, "S"))

@show time_series.u
@show time_series.w
@show time_series.T
@show time_series.S

# Coordinate arrays
xw, yw, zw = nodes(time_series.w)
xT, yT, zT = nodes(time_series.T)

times = time_series.w.times
intro = 1

n = Observable(intro)

wₙ = @lift interior(time_series.w[$n],  :, :, 1)
Tₙ = @lift interior(time_series.T[$n],  :, :, 1)
Sₙ = @lift interior(time_series.S[$n],  :, :, 1)
uₙ = @lift interior(time_series.u[$n],  :, :, 1)
vₙ = @lift interior(time_series.v[$n],  :, :, 1)

fig = Figure(resolution = (1000, 500))

axis_kwargs = (xlabel="x (m)",
               ylabel="y (m)",
               aspect = AxisAspect(grid.Lx/grid.Ly),
               limits = ((-grid.Lx/2, grid.Lx/2), (0, grid.Ly)))

ax_w  = Axis(fig[2, 1]; title = "Vertical velocity", axis_kwargs...)
ax_T  = Axis(fig[2, 3]; title = "Temperature", axis_kwargs...)
#ax_S  = Axis(fig[3, 1]; title = "Salinity", axis_kwargs...)
ax_v  = Axis(fig[3, 1]; title = "Meridional velocity", axis_kwargs...)
ax_u  = Axis(fig[3, 3]; title = "Zonal velocity", axis_kwargs...)

title = @lift @sprintf("t = %s", prettytime(times[$n]))

wlims = (-3e-6, 3e-6)
Tlims = (-1.0, 1.0)
Slims = (35, 35.005)
ulims = (-0.0005, 0.0005)


hm_w = heatmap!(ax_w, xw, yw, wₙ; colormap = :balance, colorrange = wlims)
Colorbar(fig[2, 2], hm_w; label = "m s⁻¹")

hm_T = heatmap!(ax_T, xT, yT, Tₙ; colormap = :thermal)#, colorrange = Tlims)
Colorbar(fig[2, 4], hm_T; label = "ᵒC")

#hm_S = heatmap!(ax_S, xT, yT, Sₙ; colormap = :haline)#, colorrange = Slims)
#Colorbar(fig[3, 2], hm_S; label = "g / kg")

hm_v = heatmap!(ax_v, xw, yw, vₙ; colormap = :balance, colorrange = ulims)
Colorbar(fig[3, 2], hm_v; label = "m s⁻¹")

hm_u = heatmap!(ax_u, xw, yw, uₙ; colormap = :balance, colorrange = ulims)
Colorbar(fig[3, 4], hm_u; label = "m s⁻¹")

fig[1, 1:4] = Label(fig, title, fontsize=24, tellwidth=false)

current_figure() # hide
fig

frames = intro:length(times)

@info "Making a motion picture of ocean wind mixing and convection..."

record(fig, filename * ".mp4", frames, framerate=8) do i
    n[] = i
end