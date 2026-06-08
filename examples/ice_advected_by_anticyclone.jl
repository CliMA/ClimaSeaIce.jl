# # Sea ice advected by an atmospheric anticyclone example
#
# This example simulates sea ice advected by an atmospheric anticyclone, based on
# the experiment described by [Mehlmann2021](@citet).
#
# The ice is driven by a moving anticyclonic atmospheric eddy and interacts with
# a cyclonic ocean eddy. This example demonstrates how to:
#
#   * set up a two-dimensional sea ice model with time-varying atmospheric forcing,
#   * prescribe moving atmospheric and ocean eddies,
#   * use elasto-visco-plastic rheology with Coriolis effects,
#   * save and visualize time series data.
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
# We set up a two-dimensional bounded domain:

using ClimaSeaIce
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using Oceananigans.BoundaryConditions
using Printf
using CairoMakie

arch = CPU()

L  = 512kilometers
𝓋ₒ = 0.01 # m s⁻¹ maximum ocean speed
𝓋ₐ = 30.0 # m s⁻¹ maximum atmospheric speed modifier

grid = RectilinearGrid(arch;
                       size = (128, 128),
                          x = (0, L),
                          y = (0, L),
                       halo = (7, 7),
                   topology = (Bounded, Bounded, Flat))

# ## Boundary conditions
#
# We set value boundary conditions for velocities:

u_bcs = FieldBoundaryConditions(north=ValueBoundaryCondition(0),
                                south=ValueBoundaryCondition(0))

v_bcs = FieldBoundaryConditions(west=ValueBoundaryCondition(0),
                                east=ValueBoundaryCondition(0))

# ## Ocean sea-ice stress
#
# We set up constant ocean velocities corresponding to a cyclonic eddy:

Uₒ = XFaceField(grid)
Vₒ = YFaceField(grid)

set!(Uₒ, (x, y) -> 𝓋ₒ * (2y - L) / L)
set!(Vₒ, (x, y) -> 𝓋ₒ * (L - 2x) / L)
fill_halo_regions!((Uₒ, Vₒ))

τₒ = SemiImplicitStress(uₑ=Uₒ, vₑ=Vₒ)

# ## Atmosphere-sea ice stress
#
# We set up atmospheric velocities corresponding to an anticyclonic eddy moving
# northeast. The eddy center moves with time:

@inline center(t) = 256kilometers + 51.2kilometers * t / day
@inline radius(x, y, t) = sqrt((x - center(t))^2 + (y - center(t))^2)
@inline speed(x, y, t) = 1 / 100 * exp(- radius(x, y, t) / 100kilometers)

@inline ua_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (  cosd(72) * (x - center(t)) + sind(72) * (y - center(t))) / 1000
@inline va_time(x, y, t) = - 𝓋ₐ * speed(x, y, t) * (- sind(72) * (x - center(t)) + cosd(72) * (y - center(t))) / 1000

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
# We use an elasto-visco-plastic rheology and WENO seventh order for advection
# of ice thickness and concentration. We include Coriolis effects:

dynamics = SeaIceMomentumEquation(grid;
                                  top_momentum_stress = (u=τₐu, v=τₐv),
                                  bottom_momentum_stress = τₒ,
                                  coriolis = FPlane(f=1e-4))

model = SeaIceModel(grid;
                    dynamics,
                    ice_thermodynamics = nothing, # No thermodynamics here
                    advection = WENO(order=7),
                    boundary_conditions = (u=u_bcs, v=v_bcs))

# ## Initial conditions
#
# We start with a concentration of `ℵ = 1` and an initial height field with
# perturbations around 0.3 meters:

h₀(x, y) = 0.3 + 0.005 * (sin(60 * x / 1000kilometers) + sin(30 * y / 1000kilometers))

set!(model, h = h₀)
set!(model, ℵ = 1)

# ## Running the simulation
#
# We run the model for 2 days with a 2-minute time step:

simulation = Simulation(model, Δt = 2minutes, stop_time = 2days)

# ## Time-varying wind stress
#
# We need to evolve the wind stress field in time. We set up a callback to
# update the atmospheric stress at each time step:

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
# We set up an output writer to save the ice thickness, concentration, and
# velocity fields:

h = model.ice_thickness
ℵ = model.ice_concentration
u, v = model.velocities

outputs = (; h, u, v, ℵ)

simulation.output_writers[:sea_ice] = JLD2Writer(model, outputs;
                                                 filename = "sea_ice_advected_by_anticyclone.jld2",
                                                 schedule = IterationInterval(5),
                                                 overwrite_existing = true)


# We can finally run the simulation

run!(simulation)

# ## Visualizing the results
#
# We load the saved time series and create an animation:

htimeseries = FieldTimeSeries("sea_ice_advected_by_anticyclone.jld2", "h")
utimeseries = FieldTimeSeries("sea_ice_advected_by_anticyclone.jld2", "u")
vtimeseries = FieldTimeSeries("sea_ice_advected_by_anticyclone.jld2", "v")
ℵtimeseries = FieldTimeSeries("sea_ice_advected_by_anticyclone.jld2", "ℵ")

Nt = length(htimeseries)
iter = Observable(1)

hi = @lift(interior(htimeseries[$iter], :, :, 1))
ℵi = @lift(interior(ℵtimeseries[$iter], :, :, 1))
ui = @lift(interior(utimeseries[$iter], :, :, 1))
vi = @lift(interior(vtimeseries[$iter], :, :, 1))

fig = Figure(size = (1200, 1200))
ax1 = Axis(fig[1, 1], title = "Sea ice thickness (m)")
ax2 = Axis(fig[1, 2], title = "Sea ice concentration")
ax3 = Axis(fig[2, 1], title = "Zonal velocity (m s⁻¹)")
ax4 = Axis(fig[2, 2], title = "Meridional velocity (m s⁻¹)")

heatmap!(ax1, hi, colormap = :magma, colorrange = (0.23, 0.37))
heatmap!(ax2, ℵi, colormap = Reverse(:deep), colorrange = (0.9, 1))
heatmap!(ax3, ui, colorrange = (-0.1, 0.1), colormap = :balance)
heatmap!(ax4, vi, colorrange = (-0.1, 0.1), colormap = :balance)

CairoMakie.record(fig, "sea_ice_advected_by_anticyclone.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
    @info "Rendering frame $i"
end
nothing #hide

# ![](sea_ice_advected_by_anticyclone.mp4)
#
# The animation shows how the ice is advected and deformed by the moving
# anticyclonic atmospheric eddy. Linear kinematic features (leads and ridges)
# form as the ice responds to the complex stress patterns.
