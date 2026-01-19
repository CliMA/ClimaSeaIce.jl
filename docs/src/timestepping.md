# Sea Ice Timestepping

ClimaSeaIce.jl supports two timestepping schemes for evolving the sea ice state:
Forward Euler and Split Runge-Kutta 3rd order. The default is Split Runge-Kutta 3rd order.

## Overview

Sea ice timestepping in ClimaSeaIce separates the physics into three distinct steps
that are performed sequentially within each substep:

1. **Dynamic step**: Advection of ice thickness and concentration
2. **Thermodynamic step**: Column physics (melting and freezing)
3. **Momentum step**: Evolution of ice velocities

This operator-splitting approach allows each physical process to be treated with
an appropriate numerical method.

## Forward Euler Timestepper

The `ForwardEuler` is the simplest timestepping scheme advances the state via:

```math
U^{n+1} = U^n + Δt G^n
```

where ``U^n`` is the state at timestep ``n`` and ``G^n`` is the tendency.

To use Forward Euler:

```jldoctest
using ClimaSeaIce
using Oceananigans

grid = RectilinearGrid(size=(16, 16, 1), extent=(1, 1, 1))
model = SeaIceModel(grid; timestepper = :ForwardEuler)

# output
SeaIceModel{CPU, RectilinearGrid}(time = 0 seconds, iteration = 0)
├── grid: 16×16×1 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 3×3×1 halo
├── timestepper: ForwardEulerTimeStepper
├── ice_thermodynamics: SlabThermodynamics
├── advection: Nothing
└── external_heat_fluxes:
    ├── top: Int64
    └── bottom: Int64
```

Forward Euler is first-order accurate and may require smaller time steps for stability.

## Split Runge-Kutta Timestepper (Default)

The Split Runge-Kutta 3rd order (`SplitRungeKutta3`) scheme is a 3rd-order accurate method that
performs three substeps per full time step. This is the default timestepper and
is recommended for most applications.

```jldoctest
using ClimaSeaIce
using Oceananigans

grid = RectilinearGrid(size=(16, 16, 1), extent=(1, 1, 1))
model = SeaIceModel(grid; timestepper = :SplitRungeKutta3)  # default

# output
SeaIceModel{CPU, RectilinearGrid}(time = 0 seconds, iteration = 0)
├── grid: 16×16×16 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 3×3×3 halo
├── timestepper: SplitRungeKuttaTimeStepper(3)
├── ice_thermodynamics: SlabThermodynamics
├── advection: Nothing
└── external_heat_fluxes:
    ├── top: Int64
    └── bottom: Int64
```

### Algorithm

The SplitRungeKutta3 scheme follows Oceananigans' implementation. At the beginning of each
full time step, the current state is cached via [`cache_current_fields!`](@ref ClimaSeaIce.cache_current_fields!).
Then, three substeps are performed via [`rk_substep!`](@ref ClimaSeaIce.rk_substep!), each with a different
time increment and weighting.

Each substep performs the following operations in sequence:

1. Compute advective tendencies for ice thickness ``h`` and concentration ``ℵ``
2. Update ``h`` and ``ℵ`` via [`dynamic_time_step!`](@ref ClimaSeaIce.dynamic_time_step!)
3. Apply thermodynamic fluxes (melting/freezing)
4. Update ice momentum (implicit or split-explicit)

### Advantages

- Third-order temporal accuracy
- Larger stable time steps compared to Forward Euler
- Consistent with Oceananigans' ocean model timestepping

## Choosing a Timestepper

| Timestepper | Order | Stability | Use Case |
|-------------|-------|-----------|----------|
| ForwardEuler | 1st | Requires small Δt | Simple tests, debugging |
| SplitRungeKutta3 | 3rd | Larger stable Δt | Production simulations |

For coupled ocean-sea ice simulations with ClimaOcean.jl, the SplitRungeKutta3 timestepper
is recommended as it matches the ocean model's default timestepping scheme.
