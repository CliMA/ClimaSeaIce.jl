# # Landfast sea ice
#
# This example shows sea ice grounding on a shallow shelf. Wind blows along a coast whose water
# deepens offshore, and where the ice keels are thick enough to reach the sea floor the bed arrests
# them: a stationary landfast belt forms against the coast while the pack streams freely offshore.
# Running the same configuration with and without [`LandfastBasalStress`](@ref) isolates the effect.
# This example demonstrates how to:
#
#   * give a sea-ice model a bathymetry it can ground on,
#   * switch on the Lemieux et al. (2015) basal stress,
#   * visualize the landfast belt it produces.
#
# ## Install dependencies
#
# ```julia
# using Pkg
# pkg"add Oceananigans, ClimaSeaIce, CairoMakie"
# ```

using ClimaSeaIce
using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie

# ## The physical domain
#
# The domain is a coastal strip, periodic along the coast and bounded across it. Unlike a pure
# sea-ice problem, the grid needs a vertical extent: the basal stress is driven by the water depth,
# so the model has to be told where the sea floor is.

Lx, Ly = 200kilometers, 100kilometers
Nx, Ny, Nz = 128, 64, 30

coastal_depth = 2      # m, at the coast
offshore_depth = 60    # m, at the open boundary

# The vertical grid must resolve the shelf. With cells thicker than the coastal depth the whole top
# cell would be immersed and the column depth would collapse to zero.

underlying_grid = RectilinearGrid(size = (Nx, Ny, Nz),
                                     x = (0, Lx),
                                     y = (0, Ly),
                                     z = (-offshore_depth, 0),
                                  halo = (4, 4, 4),
                              topology = (Periodic, Bounded, Bounded))

@inline sea_floor(x, y) = - (coastal_depth + (offshore_depth - coastal_depth) * y / Ly)

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(sea_floor))

# ## Where the ice grounds
#
# Keels ground where the ice is thicker than the critical thickness ``hᶜ = H ℵ / k₁``. With ``k₁ = 8``
# and three metres of fully concentrated ice, that is everywhere shallower than 24 m:

ice_thickness = 3.0
critical_thickness_parameter = 8

grounding_depth = critical_thickness_parameter * ice_thickness
grounding_line = Ly * (grounding_depth - coastal_depth) / (offshore_depth - coastal_depth)

@info @sprintf("Ice grounds in water shallower than %.0f m, i.e. within %.0f km of the coast", grounding_depth, grounding_line / 1e3)

# ## Atmospheric and oceanic forcing
#
# The wind blows along the coast with a component pushing the ice onshore, so the pack is driven
# against the shelf:

uₐ = 15.0    # m s⁻¹
Cᴰ = 1.2e-3  # atmosphere-sea ice drag coefficient
ρₐ = 1.3     # kg m⁻³

τₐ = ρₐ * Cᴰ * uₐ^2

τᵤ = XFaceField(grid)
τᵥ = YFaceField(grid)
set!(τᵤ, - τₐ)          # along the coast
set!(τᵥ, - 0.3 * τₐ)    # and towards it
Oceananigans.BoundaryConditions.fill_halo_regions!(τᵤ)
Oceananigans.BoundaryConditions.fill_halo_regions!(τᵥ)

# The ocean is at rest, so its stress opposes any ice motion:

τₒ = SemiImplicitStress()

# ## Running the two configurations
#
# The only difference between the runs is whether the sea floor is allowed to grip the ice.

function run_landfast_simulation(basal_stress; stop_time = 2days, Δt = 10minutes)
    dynamics = SeaIceMomentumEquation(grid;
                                      top_momentum_stress = (u=τᵤ, v=τᵥ),
                                      bottom_momentum_stress = τₒ,
                                      rheology = ElastoViscoPlasticRheology(),
                                      basal_stress,
                                      solver = SplitExplicitSolver(grid; substeps=150))

    model = SeaIceModel(grid; dynamics,
                        advection = WENO(order=7),
                        ice_thermodynamics = nothing)

    set!(model, h = ice_thickness)
    set!(model, ℵ = 1)

    simulation = Simulation(model; Δt, stop_time)

    speed = []
    thickness = []

    s = Field(sqrt(model.velocities.u^2 + model.velocities.v^2))

    function accumulate!(sim)
        compute!(s)
        push!(speed, deepcopy(Array(interior(s))[:, :, 1]))
        push!(thickness, deepcopy(Array(interior(sim.model.ice_thickness))[:, :, 1]))
    end

    simulation.callbacks[:save] = Callback(accumulate!, IterationInterval(6))

    run!(simulation)

    return speed, thickness
end

drifting_speed, drifting_thickness = run_landfast_simulation(nothing)
landfast_speed, landfast_thickness = run_landfast_simulation(LandfastBasalStress())

# ## Visualizing the landfast belt

x = xnodes(grid, Center()) ./ 1e3
y = ynodes(grid, Center()) ./ 1e3

Nt = length(landfast_speed)
iter = Observable(Nt)

sᵈ = @lift(100 .* drifting_speed[$iter])
sˡ = @lift(100 .* landfast_speed[$iter])
hˡ = @lift(landfast_thickness[$iter])

fig = Figure(size = (1100, 640))

ax = Axis(fig[1, 1], title = "Ice speed without basal stress (cm/s)",
          xlabel = "x (km)", ylabel = "y (km)")
heatmap!(ax, x, y, sᵈ, colormap = :speed, colorrange = (0, 30))
hlines!(ax, grounding_line / 1e3, color = :white, linestyle = :dash, linewidth = 2)

ax = Axis(fig[1, 2], title = "Ice speed with basal stress (cm/s)",
          xlabel = "x (km)", ylabel = "y (km)")
hm = heatmap!(ax, x, y, sˡ, colormap = :speed, colorrange = (0, 30))
hlines!(ax, grounding_line / 1e3, color = :white, linestyle = :dash, linewidth = 2)
Colorbar(fig[1, 3], hm)

# The cross-shore profile makes the arrest unambiguous: the belt inshore of the grounding line stops
# while the pack offshore of it keeps drifting.

ax = Axis(fig[2, 1:3], title = "Cross-shore profile of ice speed",
          xlabel = "Distance from the coast (km)", ylabel = "Speed (cm/s)")
lines!(ax, y, @lift(vec(sum($sᵈ, dims=1)) ./ Nx), linewidth = 3, label = "no basal stress")
lines!(ax, y, @lift(vec(sum($sˡ, dims=1)) ./ Nx), linewidth = 3, label = "basal stress")
vlines!(ax, grounding_line / 1e3, color = :black, linestyle = :dash, label = "grounding line")
axislegend(ax, position = :rb)
ylims!(ax, 0, 30)

CairoMakie.record(fig, "landfast_sea_ice.mp4", 1:Nt, framerate = 8) do i
    iter[] = i
end
nothing #hide

# ![](landfast_sea_ice.mp4)
#
# Without the basal stress the whole pack drifts at a uniform 25 cm/s, and the coast is simply a wall
# the ice piles against. With it, the ice inshore of the grounding line is held by its keels and stops
# dead. The step in the cross-shore profile is the landfast edge, and it falls where the bathymetry
# crosses ``H = k₁ h / ℵ`` — a line fixed by the parameters, not fitted to the result.
#
# The pack offshore of the belt is not merely left alone: it slows from 25 to about 11 cm/s. Sea ice
# is a continuum, so the grounded belt buttresses the ice beyond it through the internal stress. That
# indirect effect reaches much further than the grounding line itself, which is why a landfast
# parameterization changes an ice pack well outside the shallow water it acts in.
