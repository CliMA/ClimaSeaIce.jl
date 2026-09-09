# # A travelling top hat: incremental remapping against flux-form advection
#
# This example advects a top hat of sea ice concentration through a periodic channel at a uniform
# speed and watches what happens to the thickness recovered from the transported quantities. It
# demonstrates:
#
#   * why the sea ice model transports concentration `ℵ` and content `𝓋 = ℵh` rather than `h` itself,
#   * how the recovered thickness `h = 𝓋/ℵ` behaves at an ice edge,
#   * the difference between flux-form advection and [`IncrementalRemapping`](@ref).
#
# The exact answer is a rigid translation: the hat keeps its shape, and the thickness stays at 1.5 m
# everywhere inside it. Nothing in the problem creates or destroys ice.
#
# ## Install dependencies
#
# ```julia
# using Pkg
# pkg"add Oceananigans, ClimaSeaIce, CairoMakie"
# ```

using ClimaSeaIce
using Oceananigans
using Oceananigans.BoundaryConditions
using Printf
using CairoMakie

# ## The domain
#
# The problem is one-dimensional. We give the grid a few cells in `y` only because the model is
# two-dimensional; nothing varies across the channel. Incremental remapping reads one cell beyond its
# stencil in both directions, so the halo has to be at least two.

Nx = 200

grid = RectilinearGrid(size = (Nx, 8, 1),
                       x = (0, 1),
                       y = (0, 1),
                       z = (-1, 0),
                       halo = (5, 5, 5),
                       topology = (Periodic, Periodic, Bounded))

# ## The flow
#
# A uniform current carries the ice to the right. One revolution returns the hat to where it started.

U = 1

u = Field{Face, Center, Center}(grid)
v = Field{Center, Face, Center}(grid)

set!(u, U)
set!(v, 0)

fill_halo_regions!(u)
fill_halo_regions!(v)

# ## The initial condition
#
# A top hat of fully packed ice, 1.5 m thick, surrounded by open water. The concentration is
# discontinuous; the thickness is uniform wherever there is ice.

iced(x) = 0.3 < x < 0.6

ℵᵢ(x, y) = iced(x) ? 1.0 : 0.0
hᵢ(x, y) = iced(x) ? 1.5 : 0.0

# ## Three configurations
#
# We compare flux-form WENO with both of the model's timesteppers against incremental remapping.
# Incremental remapping transports the material swept through each face over the whole step, so it
# carries its own time integration and is a Forward-Euler scheme by construction.

configurations = [(name = "WENO5, RK3",            advection = WENO(order=5),          timestepper = :SplitRungeKutta3),
                  (name = "WENO5, Forward Euler",  advection = WENO(order=5),          timestepper = :ForwardEuler),
                  (name = "Incremental remapping", advection = IncrementalRemapping(), timestepper = :ForwardEuler),
                  (name = "UpwindBiased 1",        advection = UpwindBiased(order=1),  timestepper = :ForwardEuler)]

# ## Running
#
# The Courant number is 0.25 and we integrate for one revolution, saving a snapshot every eight steps.

courant = 0.25
Δt = courant / (U * Nx)
steps = round(Int, 1 / (U * Δt))
save_interval = 8

function run(configuration)
    model = SeaIceModel(grid; velocities = (; u, v),
                        advection = configuration.advection,
                        timestepper = configuration.timestepper,
                        ice_thermodynamics = nothing,
                        dynamics = nothing)

    set!(model, ℵ = ℵᵢ, h = hᵢ)

    ℵ = [Array(interior(model.ice_concentration))[:, 1, 1]]
    h = [Array(interior(model.ice_thickness))[:, 1, 1]]

    for n in 1:steps
        time_step!(model, Δt)

        if n % save_interval == 0
            push!(ℵ, Array(interior(model.ice_concentration))[:, 1, 1])
            push!(h, Array(interior(model.ice_thickness))[:, 1, 1])
        end
    end

    @printf("%-22s peak h = %12.4f m\n", configuration.name, maximum(maximum, h))

    return (; ℵ, h)
end

solutions = map(run, configurations)

# ## Visualizing
#
# The top panel shows the concentration, the bottom panel the recovered thickness where the
# concentration is above 10⁻³. The thickness axis is clipped at 3 m: the exact answer never leaves
# 1.5 m, so anything that runs off the top of the panel is the `h = 𝓋/ℵ` recovery failing at the
# ice edge.

x = Array(xnodes(grid, Center()))

frames = length(solutions[1].ℵ)
iter = Observable(1)

colors = [:seagreen, :crimson, :royalblue, :darkorange]

fig = Figure(size = (1000, 620))

axℵ = Axis(fig[1, 1], ylabel = "concentration ℵ", limits = ((0, 1), (-0.05, 1.15)))
axh = Axis(fig[2, 1], ylabel = "thickness h (m)", xlabel = "x", limits = ((0, 1), (0, 3)))

hlines!(axh, [1.5], color = :black, linestyle = :dash, linewidth = 1)

## The exact solution is the initial hat, translated and wrapped around the periodic domain.
time(i) = (i - 1) * save_interval * Δt
ℵexact = @lift([iced(mod(ξ - U * time($iter), 1)) ? 1.0 : 0.0 for ξ in x])
hexact = @lift([iced(mod(ξ - U * time($iter), 1)) ? 1.5 : 0.0 for ξ in x])

lines!(axℵ, x, ℵexact, color = (:black, 0.35), linewidth = 6, label = "exact")
lines!(axh, x, hexact, color = (:black, 0.35), linewidth = 6)

## The thickness is only meaningful where there is ice to carry it. Each scheme leaves a different
## trail of vanishing concentration behind the hat, and `h` is defined out to wherever that trail
## survives the model's `ℵᵐⁱⁿ` floor — so plotting `h` unmasked makes the schemes look shifted
## relative to one another when it is only their tails that differ.
ice_threshold = 1e-3

for (n, configuration) in enumerate(configurations)
    ℵn = @lift(solutions[n].ℵ[$iter])
    hn = @lift([solutions[n].ℵ[$iter][k] > ice_threshold ? solutions[n].h[$iter][k] : NaN
                for k in eachindex(x)])

    lines!(axℵ, x, ℵn, color = colors[n], linewidth = 2, label = configuration.name)
    lines!(axh, x, hn, color = colors[n], linewidth = 2)
end

axislegend(axℵ, position = :rt, framevisible = false)

title = @lift(@sprintf("Travelling top hat — t = %.2f revolutions, peak h = %.1f m",
                       U * time($iter),
                       maximum(maximum(solutions[n].h[$iter]) for n in 1:length(configurations))))

Label(fig[0, 1], title, fontsize = 18, tellwidth = false)

CairoMakie.record(fig, "one_dimensional_ice_advection.mp4", 1:frames, framerate = 12) do i
    iter[] = i
end
nothing #hide

# ![](one_dimensional_ice_advection.mp4)
#
# The concentration behaves as expected in every case: the hat translates, its edges smear — a lot for
# first-order upwinding — and it stays between zero and one.
#
# The thickness does not. Flux-form advection reconstructs `ℵ` and `𝓋 = ℵh` independently, so their
# ratio is a weighted average whose numerator and denominator come from different reconstructions.
# Where the concentration collapses at the ice edge, the denominator vanishes but the numerator does
# not vanish with it, and the recovered thickness leaves the physical range — dramatically so under
# Forward Euler, and the same mechanism is what produces spurious thick ice in a full simulation.
#
# Incremental remapping integrates the transported area and the transported content over the *same*
# swept region, so
#
# ```math
# h⁺ = \frac{Σ_k w_k ℵ̃_k h̃_k}{Σ_k w_k ℵ̃_k}
# ```
#
# is a weighted average of the reconstructed thickness with non-negative weights. It is bounded by
# construction, and here it holds 1.5 m to round-off.
