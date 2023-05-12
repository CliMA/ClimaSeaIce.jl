using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using GLMakie

# Diffusivity profile within the snow/ice
Nx = 1
Nx = 1
Nz = 10 # one snow, two ice layers

#=
grid = RectilinearGrid(size = (Nx, Ny, Nz),
                       x = (0, 100kilometers),
                       y = (0, 100kilometers),
                       z = (0, 2))
=#

grid = RectilinearGrid(size=Nz, z=(0, 2), topology=(Flat, Flat, Bounded))

vitd = VerticallyImplicitTimeDiscretization()
κ = zeros(size(grid)...)
κ .= 1e-4

closure = VerticalScalarDiffusivity(vitd, κ=κ)

ρ = 1000
Cᵖ = 4000
Qʰ = 1000 # W m⁻²
Qᵀ = Qʰ / (ρ * Cᵖ)

top_T_bc = FluxBoundaryCondition(Qᵀ)
T_bcs = FieldBoundaryConditions(top=top_T_bc)

# Note: this means we solve
#
# ∂T/∂t = ∂/∂z (κ ∂T/∂z)
model = HydrostaticFreeSurfaceModel(; grid, closure,
                                    velocities = PrescribedVelocityFields(),
                                    tracers = :T,
                                    boundary_conditions = (; T=T_bcs),
                                    buoyancy = nothing)

set!(model, T=-2)
simulation = Simulation(model, Δt=1minute, stop_iteration=3)
run!(simulation)

T = model.tracers.T
Tn = interior(T, 1, 1, :)
z = znodes(T)
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, Tn, z)
display(fig)
