
#using CairoMakie
using Oceananigans
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity # FIND submodule with CATKE
using Oceananigans.Units: minute, minutes, hour, days, kilometers
using SeawaterPolynomials.TEOS10

# This file sets up a model that resembles the Antarctic Slope Current (ASC) model in the
# 2022 paper by Si, Stewart, and Eisenman

arch = GPU()

g_Earth = 9.80665

Lx = 400kilometers
Ly = 450kilometers
Lz = 4000

Nx = 400 #64 #400
Ny = 450 #64 #450
Nz = 70 #8 #70 # TODO: modify spacing if needed, 10 m at surface, 100m at seafloor

sponge_width = 20kilometers

#
# Setting up the grid and bathymetry:
#
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

#
# We'll use Relaxation() to impose a sponge layer forcing on velocity, temperature, and salinity
#
damping_rate       = 1 / 43200 # Relaxation time scale in seconds, need to make this decrease linearly toward outermost boundary
south_mask         = GaussianMask{:y}(center=0, width=sponge_width)
north_mask         = GaussianMask{:y}(center=Ly, width=sponge_width)
south_sponge_layer = Relaxation(; rate=damping_rate, mask=south_mask)
north_sponge_layer = Relaxation(; rate=damping_rate, mask=north_mask)
sponge_layers      = (south_sponge_layer, north_sponge_layer)
# TODO: compose north_mask and south_mask together into one sponge layer, OR compose north/south sponge layers

#
# Boundary Conditions:
# (from ocean wind mixing and convection example,
# https://clima.github.io/OceananigansDocumentation/stable/generated/ocean_wind_mixing_and_convection/)
#

"""
    boundary_ramp(y, Δy)

Linear ramp from 1 to 0 between y = 0 and y = δy

For example:
```
                y=0 => boundary_ramp = 1
         0 < y < δy => boundary_ramp = (δy - y) / δy
             y > δy => boundary_ramp = 0
```
"""
δy = 10kilometers
boundary_ramp(y, δy) = min(max(0, (δy - y)/δy), 1)

u₁₀(y) = -6 * (Ly - y) / Ly # m s⁻¹, average zonal wind velocity 10 meters above the ocean at the southern boundary
v₁₀(y) =  6 * (Ly - y) / Ly # m s⁻¹, average meridional wind velocity 10 meters above the ocean at the southern boundary

cᴰ  = 2.5e-3 # dimensionless drag coefficient
ρₐ  = 1.225  # kg m⁻³, average density of air at sea-level
ρₒ  = 1026.0 # kg m⁻³, average density at the surface of the world ocean

# TODO: make this only apply at Southern boundary and decay to 0 elsewhere
Qᵘ(x, y, z) = - ρₐ / ρₒ * cᴰ * u₁₀(x) * abs(u₁₀(y)) # m² s⁻²
Qᵛ(x, y, z) = - ρₐ / ρₒ * cᴰ * v₁₀(x) * abs(v₁₀(y)) # m² s⁻²

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵛ))

# Buoyancy Equations of State - we want high order polynomials, so we'll use TEOS-10
eos = TEOS10EquationOfState() # can compare to linear EOS later (linear not recommended for polar regions)

# Coriolis Effect, using basic f-plane with precribed reference Coriolis parameter
coriolis = BetaPlane(f₀=1.3e-4, β=1e-11) #FPlane(latitude=-60)

# Diffusivities as part of closure
# TODO: make sure this works for biharmonic diffusivities as the horizontal, 
horizontal_closure = HorizontalScalarDiffusivity(ν=0.1, κ=0.1)
vertical_closure   = VerticalScalarDiffusivity(ν=3e-4, κ=1e-5)

#
# Assuming no particles or biogeochemistry
#
model = HydrostaticFreeSurfaceModel(;     grid = underlying_grid,
                            momentum_advection = WENO(),
                              tracer_advection = WENO(),
                                      buoyancy = SeawaterBuoyancy(equation_of_state=eos, gravitational_acceleration=g_Earth),
                                      coriolis = coriolis,
                                  free_surface = ImplicitFreeSurface(gravitational_acceleration=g_Earth),
                                       #forcing = (u=sponge_layers, v=sponge_layers, w=sponge_layers, T=sponge_layers, S=sponge_layers), # NamedTuple()
                                       closure = CATKEVerticalDiffusivity(),
                           #boundary_conditions = (u=u_bcs, v=v_bcs),
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

N² = 1e-8 # [s⁻²] vertical stratification, was 1e-5
M² = 1e-7 # [s⁻²] meridional temperature gradient

Δy = 100kilometers # width of the region of the front
ΔT = Δy * M²       # temperature jump associated with the front
ϵT = 1e-2 * ΔT     # noise amplitude

Tᵢ(x, y, z) = N² * z + ΔT * ramp(y, Δy) + ϵT * randn()

set!(model, T=Tᵢ)

# TODO: impose a temperature restoring at the surface that corresponds to Si/Stewart to help get baroclinic eddies
# Add initial condition to get turbulence running (should be doable on a laptop), then:
# Buoyancy gradient imposed at North/South boundary, surface temp as freezing temperature imposed

@show model

#
# Now create a simulation and run the model
#
# Full resolution is 100 sec
simulation = Simulation(model; Δt=100.0, stop_time=60days)

filename = "asc_model_60_days_Nsq_is_e-8_no_slope_hi_res_custom_beta_plane"

# Here we'll try also running a zonal average of the simulation:
u, v, w = model.velocities
avgT = Average(model.tracers.T, dims=1)
avgU = Average(u, dims=1)
avgV = Average(v, dims=1)
avgW = Average(w, dims=1)

simulation.output_writers[:zonal] = JLD2OutputWriter(model, (; T=avgT, u=avgU, v=avgV, w=avgW);
                                                     filename = filename * "_zonal_average.jld2",
                                                     schedule = IterationInterval(20),
                                                     overwrite_existing = true)


simulation.output_writers[:slices] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                     filename = filename * "_surface.jld2",
                     indices = (:, :, grid.Nz),
                     schedule = IterationInterval(20),
                     overwrite_existing = true)

run!(simulation)

@info "Simulation completed in " * prettytime(simulation.run_wall_time)
@show simulation
#=
#
# Make a figure and plot it
#

using GLMakie
using Printf

surface_filepath = filename * "_surface.jld2"
average_filepath = filename * "_zonal_average.jld2"

surface_time_series = (u = FieldTimeSeries(surface_filepath, "u"),
                       v = FieldTimeSeries(surface_filepath, "v"),
                       w = FieldTimeSeries(surface_filepath, "w"),
                       T = FieldTimeSeries(surface_filepath, "T"),
                       S = FieldTimeSeries(surface_filepath, "S"))

#@show surface_time_series.u
#@show surface_time_series.w
#@show surface_time_series.T
#@show surface_time_series.S

# Coordinate arrays
srf_xw, srf_yw, srf_zw = nodes(surface_time_series.w)
srf_xT, srf_yT, srf_zT = nodes(surface_time_series.T)

times = surface_time_series.w.times
intro = 1

n = Observable(intro)

srf_wₙ = @lift interior(surface_time_series.w[$n],  :, :, 1)
srf_Tₙ = @lift interior(surface_time_series.T[$n],  :, :, 1)
srf_Sₙ = @lift interior(surface_time_series.S[$n],  :, :, 1)
srf_uₙ = @lift interior(surface_time_series.u[$n],  :, :, 1)
srf_vₙ = @lift interior(surface_time_series.v[$n],  :, :, 1)

fig = Figure(resolution = (4000, 2000))

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

wlims = (-0.002, 0.002)
Tlims = (-0.02, 0.02)
Slims = (35, 35.005)
ulims = (-0.12, 0.12)
vlims = (-0.12, 0.12)


hm_w = heatmap!(ax_w, srf_xw, srf_yw, srf_wₙ; colormap = :balance)#, colorrange = wlims)
Colorbar(fig[2, 2], hm_w; label = "m s⁻¹")

hm_T = heatmap!(ax_T, srf_xT, srf_yT, srf_Tₙ; colormap = :thermal)#, colorrange = Tlims)
Colorbar(fig[2, 4], hm_T; label = "ᵒC")

#hm_S = heatmap!(ax_S, srf_xT, srf_yT, srf_Sₙ; colormap = :haline)#, colorrange = Slims)
#Colorbar(fig[3, 2], hm_S; label = "g / kg")

hm_v = heatmap!(ax_v, srf_xw, srf_yw, srf_vₙ; colormap = :balance)#, colorrange = vlims)
Colorbar(fig[3, 2], hm_v; label = "m s⁻¹")

hm_u = heatmap!(ax_u, srf_xw, srf_yw, srf_uₙ; colormap = :balance)#, colorrange = ulims)
Colorbar(fig[3, 4], hm_u; label = "m s⁻¹")

fig[1, 1:4] = Label(fig, title, fontsize=24, tellwidth=false)

current_figure() # hide
fig

frames = intro:length(times)

@info "Making a motion picture of ocean wind mixing and convection..."

record(fig, filename * "_surface.mp4", frames, framerate=8) do i
    n[] = i
end

#
# Now do all the above, but for zonal average time series data:
#

average_time_series = (u = FieldTimeSeries(average_filepath, "u"),
                       v = FieldTimeSeries(average_filepath, "v"),
                       w = FieldTimeSeries(average_filepath, "w"),
                       T = FieldTimeSeries(average_filepath, "T"))

#@show average_time_series.u
#@show average_time_series.w
#@show average_time_series.T
#@show average_time_series.v

# Coordinate arrays
avg_xw, avg_yw, avg_zw = nodes(average_time_series.w)
avg_xT, avg_yT, avg_zT = nodes(average_time_series.T)

times = average_time_series.w.times
intro = 1

n = Observable(intro)

avg_wₙ = @lift interior(average_time_series.w[$n],  1, :, :)
avg_Tₙ = @lift interior(average_time_series.T[$n],  1, :, :)
#avg_Sₙ = @lift interior(average_time_series.S[$n],  1, :, :)
avg_uₙ = @lift interior(average_time_series.u[$n],  1, :, :)
avg_vₙ = @lift interior(average_time_series.v[$n],  1, :, :)

fig = Figure(resolution = (4000, 2000))

axis_kwargs = (xlabel="y (m)",
               ylabel="z (m)",
               aspect = AxisAspect(2.),#AxisAspect(grid.Ly/grid.Lz),
               limits = ((0, grid.Ly), (-grid.Lz, 0)))

ax_w  = Axis(fig[2, 1]; title = "Vertical velocity", axis_kwargs...)
ax_T  = Axis(fig[2, 3]; title = "Temperature", axis_kwargs...)
#ax_S  = Axis(fig[3, 1]; title = "Salinity", axis_kwargs...)
ax_v  = Axis(fig[3, 1]; title = "Meridional velocity", axis_kwargs...)
ax_u  = Axis(fig[3, 3]; title = "Zonal velocity", axis_kwargs...)

title = @lift @sprintf("t = %s", prettytime(times[$n]))

wlims = (-0.002, 0.002)
Tlims = (-0.02, 0.02)
Slims = (35, 35.005)
ulims = (-0.12, 0.12)
vlims = (-0.12, 0.12)


hm_w = heatmap!(ax_w, avg_yw, avg_zw, avg_wₙ; colormap = :balance) #, colorrange = wlims)
Colorbar(fig[2, 2], hm_w; label = "m s⁻¹")

hm_T = heatmap!(ax_T, avg_yT, avg_zT, avg_Tₙ; colormap = :thermal) #, colorrange = Tlims)
Colorbar(fig[2, 4], hm_T; label = "ᵒC")

#hm_S = heatmap!(ax_S, srf_xT, srf_yT, srf_Sₙ; colormap = :haline)#, colorrange = Slims)
#Colorbar(fig[3, 2], hm_S; label = "g / kg")

hm_v = heatmap!(ax_v, avg_yw, avg_zw, avg_vₙ; colormap = :balance) #, colorrange = vlims)
Colorbar(fig[3, 2], hm_v; label = "m s⁻¹")

hm_u = heatmap!(ax_u, avg_yw, avg_zw, avg_uₙ; colormap = :balance) #, colorrange = ulims)
Colorbar(fig[3, 4], hm_u; label = "m s⁻¹")

fig[1, 1:4] = Label(fig, title, fontsize=24, tellwidth=false)

current_figure() # hide
fig

frames = intro:length(times)

@info "Making a motion picture of ocean wind mixing and convection..."

record(fig, filename * "_average.mp4", frames, framerate=8) do i
    n[] = i
end
=#