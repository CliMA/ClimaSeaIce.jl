
using CairoMakie
using Printf
using Oceananigans
using Oceananigans.Architectures: arch_array
using Oceananigans.Fields: ZeroField, ConstantField
using Oceananigans.AbstractOperations
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity # FIND submodule with CATKE
using Oceananigans.Units: minute, minutes, hour, days, kilometers

using SeawaterPolynomials: TEOS10EquationOfState, haline_contraction

using Statistics

using ClimaSeaIce
using ClimaSeaIce: melting_temperature
using ClimaSeaIce.HeatBoundaryConditions: RadiativeEmission, IceWaterThermalEquilibrium

include("../validation/ice_ocean_model/ice_ocean_model.jl")

# This file sets up a ocean_model that resembles the Antarctic Slope Current (ASC) ocean_model in the
# 2022 paper by Si, Stewart, and Eisenman

arch = CPU()

g_Earth = 9.80665

Lx = 400kilometers
Ly = 450kilometers
Lz = 4000

Nx = 64
Ny = 64
Nz = 8 # TODO: modify spacing if needed, 10 m at surface, 100m at seafloor

x = (-Lx/2, Lx/2)
y = (0, Ly)

topology = (Periodic, Bounded, Bounded)
halo     = (4, 4, 4)

sponge_width = 20kilometers

#
# Setting up the ocean_grid and bathymetry:
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

ice_grid = RectilinearGrid(arch,
                           size = (Nx, Ny),
                           topology = (topology[1], topology[2], Flat),
                           x = x,
                           y = y,
                           halo = halo[1:2])

underlying_ocean_grid = RectilinearGrid(arch,
                                      size = (Nx, Ny, Nz),
                                  topology = topology,
                                         x = x,
                                         y = y,
                                         z = z_faces,
                                      halo = halo)

@show ice_grid                                      
@show underlying_ocean_grid

# We want the underwater slope to provide a depth of 500 m at y = 0 and the full 4 km by y =200. It follows
# a hyperbolic tangent curve with its center at y ~= 150 at a depth of ~ -4000 + (4000 - 500)/2 = -2250 m
y₀ = 150kilometers
Δy = 25kilometers

slope_depth = 500
basin_depth = 4000

""" Varies smoothly between 0 when y << y₀ to 1 when y >> y₀ + Δy, over a region of width Δy. """
step(y, y₀, Δy) = (1 + tanh((y - y₀) / Δy)) / 2

underwater_slope(x, y) = (-basin_depth - slope_depth + (slope_depth - basin_depth) * tanh((y - y₀) / Δy)) / 2
# TODO: add 50km troughs to the slope, dependent on x and y. Use appendix C to add detail here

ocean_grid = ImmersedBoundaryGrid(underlying_ocean_grid, GridFittedBottom(underwater_slope))

@show ocean_grid

# TODO: add underwater slope at y = 0 with troughs

#
# Creating the Hydrostatic ocean_model:
#

# Forcings (ignored for now):
#=
#
# We'll use Relaxation() to impose a sponge layer forcing on velocity, temperature, and salinity
#

damping_rate       = 1 / 43200 # Relaxation time scale in seconds, need to make this decrease linearly toward outermost boundary
#south_mask         = GaussianMask{:y}(center=0, width=sponge_width)
#north_mask         = GaussianMask{:y}(center=Ly, width=sponge_width)
south_mask(x,y,z)  = y < sponge_width
north_mask(x,y,z)  = y > (Ly - sponge_width)
south_sponge_layer = Relaxation(; rate=damping_rate, mask=south_mask)
north_sponge_layer = Relaxation(; rate=damping_rate, mask=north_mask, target=0.0)#, target=LinearTarget{:z}(intercept=1.0, gradient=0))
sponge_layers      = (south_sponge_layer, north_sponge_layer)

north_sponge_layer_T = Relaxation(; rate=damping_rate, mask=north_mask, target=1.0)
sponge_layers_T      = (south_sponge_layer, north_sponge_layer_T)
=#
# TODO: compose north_mask and south_mask together into one sponge layer, OR compose north/south sponge layers

#
# Boundary Conditions:
#

# Top boundary conditions for ice/ocean:
#   - outgoing radiative fluxes emitted from surface
#   - incoming shortwave radiation starting after 40 days

ice_ocean_heat_flux      = Field{Center, Center, Nothing}(ice_grid)
top_ocean_heat_flux = Qᵀ = Field{Center, Center, Nothing}(ice_grid)
top_salt_flux       = Qˢ = Field{Center, Center, Nothing}(ice_grid)

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


@inline u₁₀(y, p) = -6 * (p.Ly - y) / p.Ly # m s⁻¹, average zonal wind velocity 10 meters above the ocean at the southern boundary
@inline v₁₀(y, p) =  6 * (p.Ly - y) / p.Ly # m s⁻¹, average meridional wind velocity 10 meters above the ocean at the southern boundary

cᴰ  = 2.5e-3 # dimensionless drag coefficient
ρₐ  = 1.225  # kg m⁻³, average density of air at sea-level
ρₒ  = 1026.0 # kg m⁻³, average density at the surface of the world ocean

# TODO: make this only apply at Southern boundary and decay to 0 elsewhere
@inline Qᵘ(x, y, t, p) = - p.ρₐ / p.ρₒ * p.cᴰ * u₁₀(y, p) * abs(u₁₀(y, p)) # m² s⁻²
@inline Qᵛ(x, y, t, p) = - p.ρₐ / p.ρₒ * p.cᴰ * v₁₀(y, p) * abs(v₁₀(y, p)) # m² s⁻²

u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ, parameters=(; cᴰ, ρₐ, ρₒ, Ly)))
v_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵛ, parameters=(; cᴰ, ρₐ, ρₒ, Ly)))

# Generate a zero-dimensional grid for a single column slab model boundary conditions:
boundary_conditions = (T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                       S = FieldBoundaryConditions(top=FluxBoundaryCondition(Qˢ)))

# Buoyancy Equations of State - we want high order polynomials, so we'll use TEOS-10
eos = TEOS10EquationOfState() # can compare to linear EOS later (linear not recommended for polar regions)

# Coriolis Effect, using basic f-plane with precribed reference Coriolis parameter
coriolis = BetaPlane(f₀=1.3e-4, β=1e-11) #FPlane(latitude=-60)

#
# Assuming no particles or biogeochemistry
#
ocean_model = HydrostaticFreeSurfaceModel(;     grid = ocean_grid,
                            momentum_advection = WENO(),
                              tracer_advection = WENO(),
                                      buoyancy = SeawaterBuoyancy(equation_of_state=eos, gravitational_acceleration=g_Earth),
                                      coriolis = coriolis,
                                  free_surface = ImplicitFreeSurface(gravitational_acceleration=g_Earth),
                                       #forcing = (u=sponge_layers, v=sponge_layers, w=sponge_layers, T=sponge_layers_T, S=sponge_layers),
                                       closure = CATKEVerticalDiffusivity(),
                           boundary_conditions = (u=u_bcs, v=v_bcs),
                                       tracers = (:T, :S, :e)
)

Nz = size(ocean_grid, 3)
So = ocean_model.tracers.S
ocean_surface_salinity = view(So, :, :, Nz)
bottom_bc = IceWaterThermalEquilibrium(ConstantField(30)) #ocean_surface_salinity)

u, v, w = ocean_model.velocities
ocean_surface_velocities = (u = view(u, :, :, Nz), #interior(u, :, :, Nz),
                            v = view(v, :, :, Nz), #interior(v, :, :, Nz),    
                            w = ZeroField())

ice_model = SlabSeaIceModel(ice_grid;
                            velocities = ocean_surface_velocities,
                            advection = nothing, #WENO(),
                            ice_consolidation_thickness = 0.05,
                            ice_salinity = 4,
                            internal_heat_flux = ConductiveFlux(conductivity=2),
                            #top_heat_flux = ConstantField(-100), # W m⁻²
                            top_heat_flux = ConstantField(0), # W m⁻²
                            top_heat_boundary_condition = PrescribedTemperature(0),
                            bottom_heat_boundary_condition = bottom_bc,
                            bottom_heat_flux = ice_ocean_heat_flux)



#
# INITIAL CONDITION FOR TEMPERATURE: a baroclinically unstable temperature distribution
#
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
M² = 1e-5 # [s⁻²] meridional temperature gradient, was 1e-7

Δy = 100kilometers # width of the region of the front
ΔT = Δy * M²       # temperature jump associated with the front
ϵT = 1e-2 * ΔT     # noise amplitude

Tᵢ(x, y, z) = N² * z + ΔT * ramp(y, Δy) + ϵT * randn()

# INITIAL CONDITION FOR ICE SALINITY:
function hᵢ(x, y)
    if sqrt(x^2 + (y-225kilometers)^2) < 100kilometers
        #return 1 + 0.1 * rand()
        return 2
    else 
        return 0
    end
end

set!(ocean_model, T=Tᵢ)
set!(ice_model, h=hᵢ)

# TODO: impose a temperature restoring at the surface that corresponds to Si/Stewart to help get baroclinic eddies
# Add initial condition to get turbulence running (should be doable on a laptop), then:
# Buoyancy gradient imposed at North/South boundary, surface temp as freezing temperature imposed

@show ice_model
@show ocean_model

#
# Now create a ocean_simulation and run the ocean_model
#

# Full resolution is 100 sec
time_step = 20minutes
stop_time = 2days

ice_simulation   = Simulation(ice_model; Δt=time_step, stop_time=stop_time, verbose=false)
ocean_simulation = Simulation(ocean_model; Δt=time_step, stop_time=stop_time, verbose=false)

filename = "asc_model_lo_res_2_days_ice"

coupled_model = IceOceanModel(ice_simulation, ocean_simulation)
coupled_simulation = Simulation(coupled_model, Δt=time_step, stop_time=stop_time)

S  = ocean_model.tracers.S
# Initial condition
S₀ = 30
T₀ = melting_temperature(ice_model.phase_transitions.liquidus, S₀) + 2.0

N²S = 1e-6
β = haline_contraction(T₀, S₀, 0, eos)
g  = ocean_model.buoyancy.model.gravitational_acceleration
by = - g * β * ∂y(S)

function progress(sim)
    h = sim.model.ice.model.ice_thickness
    S = sim.model.ocean.model.tracers.S
    T = sim.model.ocean.model.tracers.T
    u = sim.model.ocean.model.velocities.u
    msg1 = @sprintf("Iter: % 6d, time: % 12s", iteration(sim), prettytime(sim))
    msg2 = @sprintf(", max(h): %.2f", maximum(h))
    msg3 = @sprintf(", min(S): %.2f", minimum(S))
    msg4 = @sprintf(", extrema(T): (%.2f, %.2f)", minimum(T), maximum(T))
    msg5 = @sprintf(", max|∂y b|: %.2e", maximum(abs, by))
    msg6 = @sprintf(", max|u|: %.2e", maximum(abs, u))
    @info msg1 * msg2 * msg3 * msg4 * msg5 * msg6
    return nothing
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

H = ice_model.ice_thickness
T = ocean_model.tracers.T
S = ocean_model.tracers.S
u, v, w = ocean_model.velocities
η = ocean_model.free_surface.η

Ht = []
Tt = []
Ft = []
Qt = []
St = []
ut = []
vt = []
ηt = []
ζt = []
tt = []

ζ = Field(∂x(v) - ∂y(u))

function saveoutput(sim)
    compute!(ζ)
    Hn = Array(interior(H, :, :, 1))
    Fn = Array(interior(Qˢ, :, :, 1))
    Qn = Array(interior(Qᵀ, :, :, 1))
    Tn = Array(interior(T, :, :, Nz))
    Sn = Array(interior(S, :, :, Nz))
    un = Array(interior(u, :, :, Nz))
    vn = Array(interior(v, :, :, Nz))
    ηn = Array(interior(η, :, :, 1))
    ζn = Array(interior(ζ, :, :, Nz))
    push!(Ht, Hn)
    push!(Ft, Fn)
    push!(Qt, Qn)
    push!(Tt, Tn)
    push!(St, Sn)
    push!(ut, un)
    push!(vt, vn)
    push!(ηt, ηn)
    push!(ζt, ζn)
    push!(tt, time(sim))
end

coupled_simulation.callbacks[:output] = Callback(saveoutput, IterationInterval(10))

run!(coupled_simulation)

#####
##### Viz for this model with ice
#####

set_theme!(Theme(fontsize=24))

x = xnodes(ocean_grid, Center())
y = ynodes(ocean_grid, Center())

fig = Figure(resolution=(2400, 700))

axh = Axis(fig[1, 1], xlabel="x (km)", ylabel="y (km)", title="Ice thickness")
axT = Axis(fig[1, 2], xlabel="x (km)", ylabel="y (km)", title="Ocean surface temperature")
axS = Axis(fig[1, 3], xlabel="x (km)", ylabel="y (km)", title="Ocean surface salinity")
axZ = Axis(fig[1, 4], xlabel="x (km)", ylabel="y (km)", title="Ocean vorticity")

Nt = length(tt)
slider = Slider(fig[2, 1:4], range=1:Nt, startvalue=Nt)
n = slider.value

title = @lift string("Melt-driven baroclinic instability after ", prettytime(tt[$n]))
Label(fig[0, 1:3], title)

Hn = @lift Ht[$n]
Fn = @lift Ft[$n]
Tn = @lift Tt[$n]
Sn = @lift St[$n]
un = @lift ut[$n]
vn = @lift vt[$n]
ηn = @lift ηt[$n]
ζn = @lift ζt[$n]
Un = @lift mean(ut[$n], dims=1)[:]

x = x ./ 1e3
y = y ./ 1e3

Stop = view(S, :, :, Nz)
Smax = maximum(Stop)
Smin = minimum(Stop)

compute!(ζ)
ζtop = view(ζ, :, :, Nz)
ζmax = maximum(abs, ζtop)
ζlim = 2e-4 #ζmax / 2

heatmap!(axh, x, y, Hn, colorrange=(0, 1), colormap=:grays)
heatmap!(axT, x, y, Tn, colormap=:heat)
heatmap!(axS, x, y, Sn, colorrange = (29, 30), colormap=:haline)
heatmap!(axZ, x, y, ζn, colorrange=(-ζlim, ζlim), colormap=:redblue)

#heatmap!(axZ, x, y, Tn, colormap=:heat)
#heatmap!(axF, x, y, Fn)

fig #display(fig)

record(fig, filename * ".mp4", 1:Nt, framerate=48) do nn
    @info string(nn)
    n[] = nn
end


#=
ocean_simulation.output_writers[:slices] = JLD2OutputWriter(ocean_model, merge(ocean_model.velocities, ocean_model.tracers),
                                                      filename = filename * "_surface.jld2",
                                                      indices = (:, :, ocean_grid.Nz),
                                                      schedule = IterationInterval(iter_interval),
                                                      overwrite_existing = true)

ocean_simulation.output_writers[:fields] = JLD2OutputWriter(ocean_model, (; ω, s),
                                                      filename = filename * "_fields.jld2",
                                                      indices = (:, :, ocean_grid.Nz),
                                                      schedule = IterationInterval(iter_interval),
                                                      overwrite_existing = true)

run!(ocean_simulation)
@info "Simulation completed in " * prettytime(ocean_simulation.run_wall_time)
@show ocean_simulation

#
# Make a figure and plot it
#

surface_filepath = filename * "_surface.jld2"
average_filepath = filename * "_zonal_average.jld2"

surface_time_series = (u = FieldTimeSeries(surface_filepath, "u"),
                       v = FieldTimeSeries(surface_filepath, "v"),
                       w = FieldTimeSeries(surface_filepath, "w"),
                       T = FieldTimeSeries(surface_filepath, "T"),
                       S = FieldTimeSeries(surface_filepath, "S"))

#@show surface_time_series.u

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
               aspect = AxisAspect(ocean_grid.Lx/ocean_grid.Ly),
               limits = ((-ocean_grid.Lx/2, ocean_grid.Lx/2), (0, ocean_grid.Ly)))

ax_w  = Axis(fig[2, 1]; title = "Vertical velocity", axis_kwargs...)
ax_T  = Axis(fig[2, 3]; title = "Temperature", axis_kwargs...)
#ax_S  = Axis(fig[3, 1]; title = "Salinity", axis_kwargs...)
ax_v  = Axis(fig[3, 1]; title = "Meridional velocity", axis_kwargs...)
ax_u  = Axis(fig[3, 3]; title = "Zonal velocity", axis_kwargs...)

title = @lift @sprintf("t = %s", prettytime(times[$n]))

# These limits are hardcoded right now, can be changed
wlims = (-0.001, 0.001)
Tlims = (-1.1, 1.1)
Slims = (35, 35.005)
ulims = (-1.2, 1.2)
vlims = (-0.2, 0.2)

hm_w = heatmap!(ax_w, srf_xw, srf_yw, srf_wₙ; colormap = :balance, colorrange = wlims)
Colorbar(fig[2, 2], hm_w; label = "m s⁻¹")

hm_T = heatmap!(ax_T, srf_xT, srf_yT, srf_Tₙ; colormap = :thermal, colorrange = Tlims)
Colorbar(fig[2, 4], hm_T; label = "ᵒC")

#hm_S = heatmap!(ax_S, srf_xT, srf_yT, srf_Sₙ; colormap = :haline)#, colorrange = Slims)
#Colorbar(fig[3, 2], hm_S; label = "g / kg")

hm_v = heatmap!(ax_v, srf_xw, srf_yw, srf_vₙ; colormap = :balance, colorrange = vlims)
Colorbar(fig[3, 2], hm_v; label = "m s⁻¹")

hm_u = heatmap!(ax_u, srf_xw, srf_yw, srf_uₙ; colormap = :balance, colorrange = ulims)
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
# Here we'll plot vorticity and speed:
#

ω_timeseries = FieldTimeSeries(filename * "_fields.jld2", "ω")
s_timeseries = FieldTimeSeries(filename * "_fields.jld2", "s")

@show ω_timeseries
@show s_timeseries

xω, yω, zω = nodes(ω_timeseries)
xs, ys, zs = nodes(s_timeseries)
nothing # hide

set_theme!(Theme(fontsize = 24))

@info "Making a neat movie of vorticity and speed..."

fig = Figure(resolution = (4000, 2000))

axis_kwargs = (xlabel = "x",
               ylabel = "y",
               limits = ((-Lx/2, Lx/2), (0, Ly)),
               aspect = AxisAspect(1))

ax_ω = Axis(fig[2, 1]; title = "Vorticity", axis_kwargs...)
ax_s = Axis(fig[2, 3]; title = "Speed", axis_kwargs...)

n = Observable(intro)

ω = @lift interior(ω_timeseries[$n], :, :, 1)
s = @lift interior(s_timeseries[$n], :, :, 1)

hm_ω = heatmap!(ax_ω, xω, yω, ω; colormap = :balance, colorrange = (-2e-4, 2e-4))
Colorbar(fig[2, 2], hm_ω; label = "m s⁻²")
hm_s = heatmap!(ax_s, xs, ys, s; colormap = :speed, colorrange = (0, 1.5))
Colorbar(fig[2, 4], hm_s; label = "m s⁻¹")

title = @lift @sprintf("t = %s", prettytime(times[$n]))
Label(fig[1, 1:2], title, fontsize=24, tellwidth=false)

current_figure() # hide
fig

frames = intro:length(times)

@info "Making a neat animation of vorticity and speed..."

record(fig, filename * "_fields.mp4", frames, framerate=8) do i
    n[] = i
end
=#
