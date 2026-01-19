using Oceananigans.TimeSteppers: SplitRungeKuttaTimeStepper
using ClimaSeaIce.SeaIceDynamics: time_step_momentum!
using ClimaSeaIce.SeaIceThermodynamics: thermodynamic_time_step!

import Oceananigans.TimeSteppers: rk_substep!, cache_current_fields!, step_lagrangian_particles!

const RKSeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:SplitRungeKuttaTimeStepper}

step_lagrangian_particles!(model::SeaIceModel, Δτ) = nothing

"""
    cache_current_fields!(model::RKSeaIceModel)

Cache the current prognostic fields (ice thickness `h`, ice concentration `ℵ`, and any additional
tracers) into the timestepper's `Ψ⁻` storage before performing a Runge-Kutta substep.

This function is called by Oceananigans' `SplitRungeKuttaTimeStepper` at the beginning of each
full time step to store the state `Uⁿ` that is needed for computing the RK3 weighted average.

See also: [`rk_substep!`](@ref), [`dynamic_time_step!`](@ref)
"""
function cache_current_fields!(model::RKSeaIceModel)
    previous_fields = model.timestepper.Ψ⁻
    model_fields = prognostic_fields(model)
    grid = model.grid
    arch = architecture(grid)

    for name in keys(model_fields)
        Ψ⁻ = previous_fields[name]
        Ψⁿ = model_fields[name]
        parent(Ψ⁻) .= parent(Ψⁿ)
    end

    return nothing
end

"""
    rk_substep!(model::RKSeaIceModel, Δτ, callbacks)

Perform a single Runge-Kutta substep for the sea ice model, advancing the state by `Δτ`.

The substep consists of three sequential operations:

1. **Dynamic step**: Compute advective tendencies and update ice thickness `h` and
   concentration `ℵ` via [`dynamic_time_step!`](@ref).

2. **Thermodynamic step**: Apply column physics (melting/freezing) via
   `thermodynamic_time_step!`. This step is performed all at once since thermodynamics
   is local column physics.

3. **Momentum step**: Advance ice velocities using either an implicit or split-explicit
   scheme via `time_step_momentum!`.

This function is called by Oceananigans' `SplitRungeKuttaTimeStepper` for each of the
three RK3 substeps within a full time step.

Arguments
=========
- `model`: A `SeaIceModel` using `SplitRungeKuttaTimeStepper`.
- `Δτ`: The substep time increment (a fraction of the full time step `Δt`).
- `callbacks`: Callbacks to execute during the substep (currently unused).

See also: [`cache_current_fields!`](@ref), [`dynamic_time_step!`](@ref)
"""
function rk_substep!(model::RKSeaIceModel, Δτ, callbacks)

    # Compute advective tendencies and update advected tracers
    compute_tendencies!(model, Δτ)
    dynamic_time_step!(model, Δτ)

    thermodynamic_time_step!(model, model.ice_thermodynamics, Δτ)

    # This is an implicit (or split-explicit) step to advance momentum.
    time_step_momentum!(model, model.dynamics, Δτ)

    return nothing
end

"""
    dynamic_time_step!(model::RKSeaIceModel, Δt)

Update ice thickness `h` and concentration `ℵ` based on advective tendencies stored in
`model.timestepper.Gⁿ` for a Runge-Kutta substep.

Unlike the Forward Euler version, this function uses the cached previous state `Ψ⁻`
(stored by [`cache_current_fields!`](@ref)) as the base state for the update:

```math
\begin{align*}
h^{n+1} = h^n + Δt G_h^n \\\\
ℵ^{n+1} = ℵ^n + Δt G_ℵ^n
\end{align*}
```

where `hⁿ` and `ℵⁿ` are retrieved from `model.timestepper.Ψ⁻`.

The kernel `_dynamic_step_tracers!` also handles:
- Clipping negative thickness and concentration values
- Resetting concentration when thickness is zero (and vice versa)
- Ridging: when `ℵ > 1`, concentration is capped at 1 and thickness is adjusted to
  conserve ice volume

Arguments
=========
- `model`: A `SeaIceModel` using `SplitRungeKuttaTimeStepper`.
- `Δt`: The time increment for this substep.

See also: [`rk_substep!`](@ref), [`cache_current_fields!`](@ref)
"""
function dynamic_time_step!(model::RKSeaIceModel, Δt)
    grid = model.grid
    arch = architecture(grid)

    h = model.ice_thickness
    ℵ = model.ice_concentration
    hⁿ = model.timestepper.Ψ⁻.h
    ℵⁿ = model.timestepper.Ψ⁻.ℵ
    tracers = model.tracers

    Gⁿ = model.timestepper.Gⁿ

    launch!(arch, grid, :xy, _dynamic_step_tracers!, h, ℵ, hⁿ, ℵⁿ, tracers, Gⁿ, Δt)

    return nothing
end
