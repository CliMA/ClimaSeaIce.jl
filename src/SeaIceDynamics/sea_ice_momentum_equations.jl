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
                           solver = SplitExplicitSolver(150),
                           minimum_concentration = 1e-3,
                           minimum_mass = 1.0)

Constructs a `SeaIceMomentumEquation` object that controls the dynamical evolution of sea-ice momentum.
The sea-ice momentum obey the following evolution equation:

```math
    ∂u                   τₒ    τₐ
    -- + f x u = ∇ ⋅ σ + --  + -- 
    ∂t                   mᵢ    mᵢ
```
where the terms (left to right) represent (1) the time derivative of the ice velocity, (2) the coriolis force.
(3) the divergence of internal stresses, (4) the ice-ocean boundary stress, and (5) the ice-atmosphere boundary stress.

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
prognostic_fields(mom::SeaIceMomentumEquation) = prognostic_fields(mom, mom.rheology)

#####
##### Checkpointing
#####

function prognostic_state(mom::SeaIceMomentumEquation)
    pf = prognostic_fields(mom)
    return prognostic_state(pf)
end

function restore_prognostic_state!(mom::SeaIceMomentumEquation, state)
    pf = prognostic_fields(mom)
    restore_prognostic_state!(pf, state)
    return mom
end
