import Oceananigans: prognostic_state, restore_prognostic_state!
using KernelAbstractions: @kernel, @index
using RootSolvers: SecantMethod, find_zero, CompactSolution
using Oceananigans.Architectures: CPU, architecture
using Oceananigans.Fields: AbstractField, interior
using Oceananigans.Grids: ZDirection, rnode, znode
using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶜᶜᶠ, σ⁻, σⁿ
using Oceananigans.Solvers: BatchedTridiagonalSolver, solve!
using Oceananigans.TimeSteppers: SplitRungeKuttaTimeStepper
using Oceananigans.Utils: launch!

@inline on_cpu(grid) = architecture(grid) isa CPU

# Enthalpy/temperature/salinity relations and the moving two-interface vertical coordinate.
include("column_energy_relations.jl")
include("sea_ice_column_discretization.jl")

# Closures and the `ColumnEnergyThermodynamics` type, its fields, constructors, and presets.
include("column_parameterizations.jl")
include("column_energy_thermodynamics.jl")

# Per-step diagnostics, transport coefficients, and boundary-flux contributions to the solves.
include("column_diagnostics.jl")
include("column_boundary_fluxes.jl")

# Implicit energy/salinity tridiagonal solves and the moving-grid volume update.
include("column_solvers.jl")
include("column_volume_update.jl")

# Coupled time stepping and conservation diagnostics.
include("column_time_stepping.jl")
include("column_budgets.jl")
