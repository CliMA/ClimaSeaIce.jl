using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using GLMakie

# Semtner parameters
κₛ = 1e-4
κₛ = 1e-4

Nx = 16
Nx = 16
Nz = 3 # one snow, two ice layers

constant_ice_thickness = 2 # meters

ice_grid = RectilinearGrid(size = (Nx, Ny, Nz),
                           x = (0, 100kilometers),
                           y = (0, 100kilometers),
                           z = (0, constant_ice_thickness))

vitd = VerticallyImplicitTimeDiscretization()
snow_diffusivity = 1e-4 # m² s⁻¹
ice_diffusivity = [1e-4, 1e-4]
semtner_closure = VerticalScalarDiffusivity(vitd, κ=ice_diffusivity)

# Note: this means we solve
#
# ∂T/∂t = ∂/∂z (κ ∂T/∂z)
model = HydrostaticFreeSurfaceModel(; grid,
                                    velocities = PrescribedVelocityFields(),
                                    tracers = :T,
                                    buoyancy = nothing)
