using ClimaSeaIce.Rheologies
using Adapt
using KernelAbstractions: @kernel, @index
using Oceananigans: restore_prognostic_state!, prognostic_state
using Oceananigans.Architectures: architecture
using Oceananigans.Utils: launch!

struct SeaIceMomentumEquation{S, C, R, F, A, ES, FT}
    coriolis :: C
    rheology :: R
    auxiliaries :: A
    solver :: S
    free_drift :: F
    external_momentum_stresses :: ES
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
                           solver = SplitExplicitSolver(grid; substeps=150),
                           minimum_concentration = 1e-3,
                           minimum_mass = 1.0)

Constructs a `SeaIceMomentumEquation` object that controls the dynamical evolution of sea-ice momentum.
The sea-ice momentum obey the following evolution equation:

```math
\\frac{∂\\boldsymbol{u}}{∂t} + \\boldsymbol{f} × \\boldsymbol{u} = \\frac{\\boldsymbol{\\nabla} \\cdot \\boldsymbol{\\sigma}}{mᵢ} + \\frac{\\boldsymbol{\\tau}ₒ}{mᵢ} + \\frac{\\boldsymbol{\\tau}ₐ}{mᵢ}
```
where ``∂\\boldsymbol{u}/∂t`` is the time derivative of the ice velocity, ``\\boldsymbol{f}`` is the
Coriolis parameter, ``\\boldsymbol{\\nabla} \\cdot \\boldsymbol{\\sigma} / mᵢ`` is the divergence of internal
stresses, ``\\boldsymbol{\\tau}ₒ/mᵢ`` is the ice-ocean boundary stress, ``\\boldsymbol{\\tau}ₐ/mᵢ`` is the
ice-atmosphere boundary stress, and ``mᵢ = ρᵢ h ℵ`` is the ice mass per unit area.

Arguments
=========

- `grid`: The computational grid.

Keyword Arguments
=================

- `coriolis`: Parameters for the background rotation rate of the model.
- `rheology`: The sea ice rheology model, default is `ElastoViscoPlasticRheology(eltype(grid))`.
- `free_drift`: The free drift velocities used when nonzero sea ice mass or concentration are below
                the dynamical momentum thresholds. Default is `nothing`.
- `solver`: The momentum solver to be used.
- `minimum_concentration`: The minimum sea ice concentration above which the sea ice velocity is dynamically calculated; below this threshold nonzero sea ice moves with free drift, and roundoff-level concentration cells are set to zero. Default is `1e-3`.
- `minimum_mass`: The minimum sea ice mass per area above which the sea ice velocity is dynamically calculated; below this threshold nonzero sea ice moves with free drift, and roundoff-level mass cells are set to zero. Default is `1.0 kg/m²`.
"""
function SeaIceMomentumEquation(grid;
                                coriolis = nothing,
                                rheology = ElastoViscoPlasticRheology(eltype(grid)),
                                top_momentum_stress    = nothing,
                                bottom_momentum_stress = nothing,
                                free_drift = nothing,
                                solver = SplitExplicitSolver(grid; substeps=150),
                                minimum_concentration = 1e-3,
                                minimum_mass = 1.0)

    auxiliaries = Auxiliaries(rheology, grid)
    external_momentum_stresses = (top = materialize_stress(top_momentum_stress, grid),
                                  bottom = materialize_stress(bottom_momentum_stress, grid))

    FT = eltype(grid)

    return SeaIceMomentumEquation(coriolis,
                                  rheology,
                                  auxiliaries,
                                  solver,
                                  free_drift,
                                  external_momentum_stresses,
                                  convert(FT, minimum_concentration),
                                  convert(FT, minimum_mass))
end

fields(mom::SeaIceMomentumEquation) = mom.auxiliaries.fields
prognostic_fields(model, mom::SeaIceMomentumEquation) = merge(model.velocities, prognostic_fields(mom, mom.rheology))

# Fallback: keep the same grid and dynamics
maybe_extended_grid(mom, grid) = grid
materialize_solver(mom, velocity_grid) = mom

@kernel function _reconcile_external_stress_coefficients!(grid, clock, fields, top_stress, bottom_stress)
    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)

    compute_implicit_stress_coefficients!(i, j, kᴺ, grid, top_stress, clock, fields)
    compute_implicit_stress_coefficients!(i, j, kᴺ, grid, bottom_stress, clock, fields)
end

reconcile_dynamics!(model, ::Nothing) = nothing
reconcile_dynamics!(model, mom::SeaIceMomentumEquation) = nothing

function Base.show(io::IO, sime::SeaIceMomentumEquation)

    aux_fields = keys(sime.auxiliaries.fields)

    print(io, "SeaIceMomentumEquation", '\n')
    print(io, "├── coriolis: ", summary(sime.coriolis), '\n')
    print(io, "├── rheology: ", summary(sime.rheology), '\n')
    print(io, "├── auxiliaries: ", join(aux_fields, ", "), '\n')
    print(io, "├── solver: ", summary(sime.solver), '\n')
    print(io, "├── free_drift: ", sime.free_drift, '\n')
    print(io, "├── external_momentum_stresses: ", keys(sime.external_momentum_stresses), '\n')
    print(io, "├── minimum_concentration: ", sime.minimum_concentration, '\n')
    print(io, "└── minimum_mass: ", sime.minimum_mass)
end

#####
##### Checkpointing
#####

function Oceananigans.prognostic_state(mom::SeaIceMomentumEquation)
    return (; fields = prognostic_state(fields(mom)))
end

function Oceananigans.restore_prognostic_state!(mom::SeaIceMomentumEquation, state)
    restore_prognostic_state!(fields(mom), state.fields)
    return mom
end
