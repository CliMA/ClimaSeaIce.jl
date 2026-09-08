using ..Rheologies

struct SeaIceMomentumEquation{S, C, R, F, A, ES, B, H, FT}
    coriolis :: C
    rheology :: R
    auxiliaries :: A
    solver :: S
    free_drift :: F
    external_momentum_stresses :: ES
    basal_stress :: B
    free_surface :: H
    minimum_concentration :: FT
    minimum_mass :: FT
end

""" a simple explicit solver """
struct ExplicitSolver end

"""
    SeaIceMomentumEquation(grid;
                           coriolis = nothing,
                           rheology = ElastoViscoPlasticRheology(eltype(grid)),
                           top_momentum_stress    = nothing,
                           bottom_momentum_stress = nothing,
                           free_drift = nothing,
                           basal_stress = nothing,
                           ocean_surface_height = ZeroField(),
                           gravitational_acceleration = Oceananigans.defaults.gravitational_acceleration,
                           solver = SplitExplicitSolver(grid; substeps=150),
                           minimum_concentration = 1e-3,
                           minimum_mass = 1.0)

Constructs a `SeaIceMomentumEquation` object that controls the dynamical evolution of sea-ice momentum.
The sea-ice momentum obey the following evolution equation:

```math
\\frac{∂\\boldsymbol{u}}{∂t} + \\boldsymbol{f} × \\boldsymbol{u} = \\frac{\\boldsymbol{\\nabla} \\cdot \\boldsymbol{\\sigma}}{mᵢ} + \\frac{\\boldsymbol{\\tau}ₒ}{mᵢ} + \\frac{\\boldsymbol{\\tau}ₐ}{mᵢ} - g \\boldsymbol{\\nabla} η
```
where ``∂\\boldsymbol{u}/∂t`` is the time derivative of the ice velocity, ``\\boldsymbol{f}`` is the
Coriolis parameter, ``\\boldsymbol{\\nabla} \\cdot \\boldsymbol{\\sigma} / mᵢ`` is the divergence of internal
stresses, ``\\boldsymbol{\\tau}ₒ/mᵢ`` is the ice-ocean boundary stress, ``\\boldsymbol{\\tau}ₐ/mᵢ`` is the
ice-atmosphere boundary stress, ``g \\boldsymbol{\\nabla} η`` is the acceleration due to the tilt of the
ocean surface, and ``mᵢ = ρᵢ h ℵ`` is the ice mass per unit area.

Arguments
=========

- `grid`: The computational grid.

Keyword Arguments
=================

- `coriolis`: Parameters for the background rotation rate of the model.
- `rheology`: The sea-ice rheology model. Default:
              `ElastoViscoPlasticRheology(eltype(grid))`.
- `top_momentum_stress`: Atmosphere-to-ice momentum stress, or an object that can
                         be materialized into one. Default: `nothing`.
- `bottom_momentum_stress`: Ocean-to-ice momentum stress, or an object that can
                            be materialized into one. Default: `nothing`.
- `free_drift`: The free drift velocities used when nonzero sea ice mass or concentration are below
                the dynamical momentum thresholds. Default is `nothing`.
- `basal_stress`: Stress exerted by the sea floor on grounded keels, arresting landfast ice over
                  shallow bathymetry. Default: `nothing`. See [`LandfastBasalStress`](@ref).
- `ocean_surface_height`: Surface height ``η`` [m] of the underlying ocean, whose slope accelerates the
                          ice down the ocean's dynamic topography. Default: `ZeroField()`.
- `gravitational_acceleration`: ``g`` (m s⁻²). Default: `Oceananigans.defaults.gravitational_acceleration`.
- `solver`: Momentum solver used to advance the velocity field. Default:
            `SplitExplicitSolver(grid; substeps = 150)`.
- `minimum_concentration`: Minimum sea-ice concentration above which the velocity
                           is evolved dynamically. Below this threshold, nonzero
                           ice moves with free drift and roundoff-level
                           concentration cells are set to zero. Default: `1e-3`.
- `minimum_mass`: Minimum sea-ice mass per area above which the velocity is
                  evolved dynamically. Below this threshold, nonzero ice moves
                  with free drift and roundoff-level mass cells are set to zero.
                  Default: `1.0 kg/m²`.
"""
function SeaIceMomentumEquation(grid;
                                coriolis = nothing,
                                rheology = ElastoViscoPlasticRheology(eltype(grid)),
                                top_momentum_stress    = nothing,
                                bottom_momentum_stress = nothing,
                                free_drift = nothing,
                                basal_stress = nothing,
                                ocean_surface_height = ZeroField(),
                                gravitational_acceleration = Oceananigans.defaults.gravitational_acceleration,
                                solver = SplitExplicitSolver(grid; substeps=150),
                                minimum_concentration = 1e-3,
                                minimum_mass = 1.0)

    auxiliaries = Auxiliaries(rheology, grid)
    external_momentum_stresses = (top = materialize_stress(top_momentum_stress, grid),
                                  bottom = materialize_stress(bottom_momentum_stress, grid))

    # Keep the free drift pointing at the same (materialized) stress fields as the external stresses.
    free_drift = materialize_free_drift(free_drift, external_momentum_stresses.top, external_momentum_stresses.bottom)
    basal_stress = materialize_basal_stress(basal_stress, grid)
    free_surface = materialize_free_surface(ocean_surface_height, gravitational_acceleration, grid)

    FT = eltype(grid)

    return SeaIceMomentumEquation(coriolis,
                                  rheology,
                                  auxiliaries,
                                  solver,
                                  free_drift,
                                  external_momentum_stresses,
                                  basal_stress,
                                  free_surface,
                                  convert(FT, minimum_concentration),
                                  convert(FT, minimum_mass))
end

Oceananigans.fields(mom::SeaIceMomentumEquation) = mom.auxiliaries.fields
Oceananigans.prognostic_fields(model, mom::SeaIceMomentumEquation) = merge(model.velocities, prognostic_fields(mom, mom.rheology))

# Fallback: keep the same grid and dynamics
maybe_extended_grid(mom, grid) = grid
materialize_solver(mom, velocity_grid) = mom

function Base.show(io::IO, sime::SeaIceMomentumEquation)

    aux_fields = keys(sime.auxiliaries.fields)

    print(io, "SeaIceMomentumEquation", '\n')
    print(io, "├── coriolis: ", summary(sime.coriolis), '\n')
    print(io, "├── rheology: ", summary(sime.rheology), '\n')
    print(io, "├── auxiliaries: ", join(aux_fields, ", "), '\n')
    print(io, "├── basal_stress: ", summary(sime.basal_stress), '\n')
    print(io, "├── ocean_surface_height: ", summary(sime.free_surface.η), '\n')
    print(io, "├── solver: ", summary(sime.solver), '\n')
    print(io, "├── free_drift: ", sime.free_drift, '\n')
    print(io, "├── external_momentum_stresses: ", keys(sime.external_momentum_stresses), '\n')
    print(io, "├── minimum_concentration: ", sime.minimum_concentration, '\n')
    print(io, "└── minimum_mass: ", sime.minimum_mass)
end

#####
##### Checkpointing
#####

Oceananigans.prognostic_state(mom::SeaIceMomentumEquation) =
    (; fields = prognostic_state(fields(mom)))

function Oceananigans.restore_prognostic_state!(mom::SeaIceMomentumEquation, state)
    restore_prognostic_state!(fields(mom), state.fields)
    return mom
end

Oceananigans.restore_prognostic_state!(mom::SeaIceMomentumEquation, ::Nothing) = nothing
