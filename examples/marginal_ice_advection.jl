# # Marginal Ice advection
#
# This example simulates an (initially) square block of ice moving with a uniform
# ocean current and deformed by a passing atmospheric cyclone overhead. The
# scenario mimics a marginal ice zone, where a discrete patch of ice is
# surrounded by open water, and demonstrates how to:
#
#   * set up a marginal ice zone with partial ice cover,
#   * combine a uniform ocean current with time-varying atmospheric forcing,
#   * model a translating cyclonic eddy passing over the ice,
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
# We set up a two-dimensional periodic channel:

using ClimaSeaIce
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Printf
using CairoMakie

Lx = 512kilometers
Ly = 256kilometers
Nx = 256
Ny = 128

arch = CPU()

grid = RectilinearGrid(arch; size = (Nx, Ny),
                                x = (-Lx/2, Lx/2),
                                y = (0, Ly),
                             halo = (7, 7),
                         topology = (Periodic, Bounded, Flat))

# ## Boundary conditions
#
# We close the northern and southern walls of the channel with no-slip
# conditions on the zonal ice velocity:

u_bcs = FieldBoundaryConditions(north = ValueBoundaryCondition(0),
                                south = ValueBoundaryCondition(0))

# ## Ocean sea-ice stress
#
# A uniform eastward ocean current provides the background advection that
# carries the marginal ice block across the channel:

𝓋ₒ = 0.2 # m s⁻¹ uniform ocean speed

Uₒ = XFaceField(grid)
set!(Uₒ, 𝓋ₒ)
fill_halo_regions!(Uₒ)

τₒ = SemiImplicitStress(uₑ=Uₒ)

# ## Atmosphere-sea ice stress
#
# The atmosphere consists in a cyclonic eddy translating eastward across the
# domain. The translation speed is chosen larger than the ocean current so
# that the cyclone overtakes the ice block during the simulation:

𝓋ₐ = 25.0 # m s⁻¹ maximum atmospheric speed modifier
y_cyclone = Ly / 2

@inline x_cyclone(t) = -Lx/2 + 40kilometers + 128kilometers * t / day
@inline cyclone_radius(x, y, t) = sqrt((x - x_cyclone(t))^2 + (y - y_cyclone)^2)
@inline speed(x, y, t) = 1 / 100 * exp(- cyclone_radius(x, y, t) / 80kilometers)

@inline ua_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (  cosd(72) * (x - x_cyclone(t)) + sind(72) * (y - y_cyclone)) / 1000
@inline va_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (- sind(72) * (x - x_cyclone(t)) + cosd(72) * (y - y_cyclone)) / 1000

# Initialize the stress at time t = 0:

Uₐ = XFaceField(grid)
Vₐ = YFaceField(grid)

set!(Uₐ, (x, y) -> ua_time(x, y, 0))
set!(Vₐ, (x, y) -> va_time(x, y, 0))

fill_halo_regions!((Uₐ, Vₐ))

τₐu = Field(- Uₐ * sqrt(Uₐ^2 + Vₐ^2) * 1.3 * 1.2e-3)
τₐv = Field(- Vₐ * sqrt(Uₐ^2 + Vₐ^2) * 1.3 * 1.2e-3)

# ## Model configuration
#
# We use an elasto-visco-plastic rheology with split-explicit time stepping
# and WENO seventh order advection for ice thickness and concentration. We
# include Coriolis effects appropriate to mid-latitudes:

dynamics = SeaIceMomentumEquation(grid;
                                  top_momentum_stress = (u=τₐu, v=τₐv),
                                  bottom_momentum_stress = τₒ,
                                  coriolis = FPlane(f=1e-4),
                                  rheology = ElastoViscoPlasticRheology(),
                                  solver = SplitExplicitSolver(substeps=150))

model = SeaIceModel(grid;
                    dynamics,
                    advection = WENO(order=7),
                    timestepper = :ForwardEuler,
                    boundary_conditions = (; u=u_bcs),
                    ice_thermodynamics = nothing)

# ## Initial conditions
#
# We initialize a square block of ice in the middle of the channel with
# smoothed edges to limit Gibbs oscillations at the ice edge. The same mask
# is used for both ice thickness and concentration:

block_radius = 64kilometers
edge_width   = 4kilometers

@inline ice_block(x, y) = 0.25 * (1 - tanh((abs(x)          - block_radius) / edge_width)) *
                                 (1 - tanh((abs(y - Ly/2)   - block_radius) / edge_width))

set!(model, h = ice_block, ℵ = ice_block)

# ## Running the simulation
#
# We run the model for 3 days with a 2-minute time step:

simulation = Simulation(model, Δt = 5minutes, stop_time = 5days)

# ## Time-varying wind stress
#
# A callback updates the atmospheric stress at each time step as the cyclone
# translates eastward across the channel:

function compute_wind_stress(sim)
    time = sim.model.clock.time
    @inline ua(x, y) = ua_time(x, y, time)
    @inline va(x, y) = va_time(x, y, time)
    set!(Uₐ, ua)
    set!(Vₐ, va)

    fill_halo_regions!((Uₐ, Vₐ))

    compute!(τₐu)
    compute!(τₐv)

    return nothing
end

simulation.callbacks[:top_stress] = Callback(compute_wind_stress, IterationInterval(1))

# ## Output writer
#
# We save the ice thickness, concentration, and velocity fields:

h = model.ice_thickness
ℵ = model.ice_concentration
u, v = model.velocities

outputs = (; h, u, v, ℵ)

simulation.output_writers[:sea_ice] = JLD2Writer(model, outputs;
                                                 filename = "marginal_ice_advection.jld2",
                                                 schedule = IterationInterval(10),
                                                 overwrite_existing = true)

# We can finally run the simulation

run!(simulation)

# ## Visualizing the results
#
# We load the saved time series and create an animation:

htimeseries = FieldTimeSeries("marginal_ice_advection.jld2", "h")
utimeseries = FieldTimeSeries("marginal_ice_advection.jld2", "u")
vtimeseries = FieldTimeSeries("marginal_ice_advection.jld2", "v")
ℵtimeseries = FieldTimeSeries("marginal_ice_advection.jld2", "ℵ")

Nt = length(htimeseries)
iter = Observable(1)

hi = @lift(interior(htimeseries[$iter], :, :, 1))
ℵi = @lift(interior(ℵtimeseries[$iter], :, :, 1))
ui = @lift(interior(utimeseries[$iter], :, :, 1))
vi = @lift(interior(vtimeseries[$iter], :, :, 1))

fig = Figure(size = (1200, 800))
ax1 = Axis(fig[1, 1], title = "Sea ice thickness (m)")
ax2 = Axis(fig[1, 2], title = "Sea ice concentration")
ax3 = Axis(fig[2, 1], title = "Zonal velocity (m s⁻¹)")
ax4 = Axis(fig[2, 2], title = "Meridional velocity (m s⁻¹)")

heatmap!(ax1, hi, colormap = :magma,          colorrange = (0.0, 1.5))
heatmap!(ax2, ℵi, colormap = Reverse(:deep),  colorrange = (0.0, 1.0))
heatmap!(ax3, ui, colormap = :balance,        colorrange = (-0.2, 0.2))
heatmap!(ax4, vi, colormap = :balance,        colorrange = (-0.2, 0.2))

CairoMakie.record(fig, "marginal_ice_advection.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
    @info "Rendering frame $i"
end
nothing #hide

# ![](marginal_ice_advection.mp4)
#
# The animation shows the ice block being advected eastward by the uniform
# ocean current while being progressively deformed as the cyclone passes
# overhead. Ridges and leads form along the ice edges as the cyclone's
# swirling winds compress and shear the marginal ice cover.
