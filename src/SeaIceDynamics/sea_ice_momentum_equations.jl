using ClimaSeaIce.Rheologies
using Adapt

import Oceananigans: prognostic_state, restore_prognostic_state!

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
\frac{∂𝐮}{∂t} + 𝐟 × 𝐮 = 𝛁 ⋅ 𝛔 + \frac{τₒ}{mᵢ} + + \frac{τₐ}{mᵢ}
```
where ``∂𝐮/∂t`` is the time derivative of the ice velocity, ``𝐟`` is the coriolis force,
``𝛁 ⋅ 𝛔`` is the divergence of internal stresses, ``τₒ/mᵢ`` is the ice-ocean boundary stress,
and ``τₐ/mᵢ`` is the ice-atmosphere boundary stress.

Arguments
=========

- `grid`: The computational grid.

Keyword Arguments
=================

- `coriolis`: Parameters for the background rotation rate of the model.
- `rheology`: The sea ice rheology model, default is `ElastoViscoPlasticRheology(eltype(grid))`.
- `free_drift`: The free drift velocities used to limit sea ice momentum when the mass or the concentration are
                below a certain threshold. Default is `nothing` (indicating that the free drift velocities are zero).
- `solver`: The momentum solver to be used.
- `minimum_concentration`: The minimum sea ice concentration above which the sea ice velocity is dynamically calculated, default is `1e-3`.
- `minimum_mass`: The minimum sea ice mass per area above which the sea ice velocity is dynamically calculated, default is `1.0 kg/m²`.
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
    external_momentum_stresses = (top = top_momentum_stress,
                                  bottom = bottom_momentum_stress)

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

function prognostic_state(mom::SeaIceMomentumEquation)
    return (; fields = prognostic_state(fields(mom)))
end

function restore_prognostic_state!(mom::SeaIceMomentumEquation, state)
    restore_prognostic_state!(fields(mom), state.fields)
    return mom
end
