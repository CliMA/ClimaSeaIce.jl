# Column Energy Thermodynamics

```@meta
CurrentModule = ClimaSeaIce.SeaIceThermodynamics
```

The column energy thermodynamics model represents a vertical sea-ice column with
prognostic internal energy and, optionally, prognostic bulk salinity. It is
designed to cover two useful subcases with the same field and solver machinery:

- fixed-salinity enthalpy thermodynamics, constructed with
  [`prescribed_salinity_enthalpy_thermodynamics`](@ref);
- evolving-salinity mushy thermodynamics, constructed with
  [`evolving_salinity_mushy_thermodynamics`](@ref).

The derivation follows the local theory note
`A_hierarchy_of_thermodynamic_sea_ice_models.pdf`: volume-averaged sea ice is
treated as a mixture of solid ice and liquid brine in microscopic thermal
equilibrium, with vanishing solid-ice salinity.

## Mixture Thermodynamics

The liquidus is linear,

```math
T_m(S) = T_0 - m S ,
```

so local thermal equilibrium implies a brine salinity

```math
S_b = \frac{T_0 - T}{m}.
```

With zero solid-ice salinity, bulk salinity is carried by the brine and the
liquid fraction is

```math
\phi = \frac{m S}{T_0 - T}.
```

The implemented energy relation is

```math
E(T, S) =
\rho_i c_i (T - T_0)
- (\rho_l c_l - \rho_i c_i) m S
- \rho_l L_0 \frac{m S}{T - T_0}.
```

For prescribed internal energy and bulk salinity this becomes a quadratic
equation for ``T - T_0``. [`QuadraticLiquidusEnergyRelation`](@ref) evaluates
both directions, plus the derivatives ``\partial T / \partial E`` and
``\partial T / \partial S`` used in the semi-implicit solve.

```@example column_energy_relation
using ClimaSeaIce.SeaIceThermodynamics:
    QuadraticLiquidusEnergyRelation,
    internal_energy,
    temperature,
    liquid_fraction,
    brine_salinity

relation = QuadraticLiquidusEnergyRelation(Float64)
E = internal_energy(relation, -10.0, 5.0)
recovered_temperature = temperature(relation, E, 5.0)
phi = liquid_fraction(relation, E, 5.0)
brine_S = brine_salinity(relation, E, 5.0)

abs(recovered_temperature + 10) < 1e-12 && 0 <= phi <= 1 && isfinite(brine_S)
```

## Transport Closures

The column solves

```math
\partial_t E = \partial_z Q_E + \partial_z I,
\qquad
\partial_t S = \partial_z J_S.
```

The currently implemented energy flux closures are

```math
Q_E = k \partial_z T,
\qquad
Q_E = \kappa_E \partial_z E,
\qquad
Q_E = k \partial_z T + \kappa_E \partial_z E,
```

represented by [`ConductiveTemperatureTransport`](@ref),
[`DiffusiveEnergyTransport`](@ref), and
[`ConductiveAndDiffusiveEnergyTransport`](@ref). Bulk-salinity transport is
either disabled with [`NoSalinityTransport`](@ref), or stepped with closed
boundary scalar diffusion via [`BulkSalinityDiffusion`](@ref).
[`BrineSalinityDiffusion`](@ref) is reserved as an explicit future marker for a
brine-salinity transport equation and is not used by the scalar bulk-salinity
step.

Boundary energy fluxes are configured with [`ColumnBoundaryConditions`](@ref).
[`InsulatingBoundary`](@ref) imposes zero flux. [`PrescribedEnergyFlux`](@ref)
imposes a face flux that is positive in the increasing vertical-coordinate
direction. [`MeltingLimitedSurfaceFlux`](@ref) imposes a top surface flux but
caps the applied fixed-grid energy at the complete-melt threshold of the top
cell, returning the excess as a Stefan residual for surface melt.
[`ExponentialShortwaveAbsorption`](@ref) adds a Beer-law shortwave flux ``I``
with prescribed transmitted flux at the top face and an e-folding attenuation
scale.

## Semi-Implicit Step

Conductive temperature transport is linearized using the thermodynamic
derivatives:

```math
\partial_z (k \partial_z T)
\approx
\partial_z \left[
\left(k \frac{\partial T}{\partial E} + \kappa_E\right) \partial_z E^{n+1}
+ k \frac{\partial T}{\partial S} \partial_z S^n
\right].
```

This yields a backward-Euler tridiagonal system for ``E^{n+1}``. The scalar
bulk-salinity diffusion step uses the same fixed-grid tridiagonal machinery when
bulk salinity is prognostic. The CPU path uses an allocation-free Thomas sweep
over the same coefficient fields; non-CPU architectures use the Oceananigans
batched tridiagonal solver path.

```@example column_energy_step
using Oceananigans
using Oceananigans.Fields: set!
using ClimaSeaIce.SeaIceThermodynamics:
    ColumnBoundaryConditions,
    ConductiveTemperatureTransport,
    InsulatingBoundary,
    prescribed_salinity_enthalpy_thermodynamics,
    column_energy_budget,
    column_energy_time_step!,
    column_integrated_energy

grid = RectilinearGrid(size = 8,
                       z = (0, 1),
                       topology = (Flat, Flat, Bounded))

thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
    salinity_profile = 0.0,
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
    boundary_conditions = ColumnBoundaryConditions(top = InsulatingBoundary(),
                                                   bottom = InsulatingBoundary()))

set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> -12 + 4z)
initial_energy = column_integrated_energy(thermodynamics)
column_energy_time_step!(thermodynamics, 5e3)
budget = column_energy_budget(thermodynamics, initial_energy, 5e3)

budget.relative_residual < 1e-12
```

Shortwave absorption enters the same budget through the face-flux difference
``I_\mathrm{top} - I_\mathrm{bottom}``.

```@example column_energy_step
using ClimaSeaIce.SeaIceThermodynamics:
    ExponentialShortwaveAbsorption,
    compute_column_shortwave_flux!

shortwave = ExponentialShortwaveAbsorption(surface_transmission = 3.0,
                                           attenuation_scale = 0.25)

shortwave_thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
    salinity_profile = 0.0,
    energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
    shortwave_absorption = shortwave,
    boundary_conditions = ColumnBoundaryConditions(top = InsulatingBoundary(),
                                                   bottom = InsulatingBoundary()))

set!(shortwave_thermodynamics; bulk_salinity = 0.0, temperature = -10.0)
initial_energy = column_integrated_energy(shortwave_thermodynamics)
column_energy_time_step!(shortwave_thermodynamics, 100.0)
budget = column_energy_budget(shortwave_thermodynamics, initial_energy, 100.0)

budget.relative_residual < 1e-11
```

## Evolving Salinity

The evolving-salinity preset adds `bulk_salinity` to the prognostic fields. With
[`NoSalinityTransport`](@ref), it reduces to the fixed-salinity model for the
same initial salinity field. With [`BulkSalinityDiffusion`](@ref), closed
boundaries conserve column-integrated salinity while reducing salinity variance.

```@example column_salinity_step
using Oceananigans
using Oceananigans.Fields: set!, interior
using ClimaSeaIce.SeaIceThermodynamics:
    BulkSalinityDiffusion,
    ColumnBoundaryConditions,
    ConductiveTemperatureTransport,
    InsulatingBoundary,
    evolving_salinity_mushy_thermodynamics,
    column_integrated_salinity,
    column_salt_budget,
    column_salinity_time_step!

grid = RectilinearGrid(size = 8,
                       z = (0, 1),
                       topology = (Flat, Flat, Bounded))

thermodynamics = evolving_salinity_mushy_thermodynamics(grid;
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
    salinity_transport = BulkSalinityDiffusion(diffusivity = 1e-4),
    boundary_conditions = ColumnBoundaryConditions(top = InsulatingBoundary(),
                                                   bottom = InsulatingBoundary()))

set!(thermodynamics;
     bulk_salinity = z -> 5 + sin(2pi * z),
     temperature = -10.0)

salinity_variance(thermodynamics) = begin
    S = vec(Array(interior(thermodynamics.fields.bulk_salinity)))
    S_mean = sum(S) / length(S)
    sum(abs2, S .- S_mean) / length(S)
end

initial_salt = column_integrated_salinity(thermodynamics)
initial_variance = salinity_variance(thermodynamics)
column_salinity_time_step!(thermodynamics, 100)
budget = column_salt_budget(thermodynamics, initial_salt, 100)

budget.relative_residual < 1e-12 &&
    salinity_variance(thermodynamics) <= initial_variance
```

## Fixed and Evolving Presets

The fixed-salinity preset is a strict subcase of the evolving-salinity preset
when salinity transport is disabled. A side-by-side run with the same initial
temperature and salinity profiles reports the resulting equivalence diagnostics.

```@example column_equivalence
using Oceananigans
using Oceananigans.Fields: set!, interior
using ClimaSeaIce.SeaIceThermodynamics:
    ColumnBoundaryConditions,
    ConductiveTemperatureTransport,
    InsulatingBoundary,
    NoSalinityTransport,
    column_energy_time_step!,
    evolving_salinity_mushy_thermodynamics,
    prescribed_salinity_enthalpy_thermodynamics

grid = RectilinearGrid(size = 8,
                       z = (0, 1),
                       topology = (Flat, Flat, Bounded))

boundary_conditions = ColumnBoundaryConditions(top = InsulatingBoundary(),
                                               bottom = InsulatingBoundary())
energy_transport = ConductiveTemperatureTransport(conductivity = 2.0)

fixed = prescribed_salinity_enthalpy_thermodynamics(grid;
    salinity_profile = 0.0,
    energy_transport,
    boundary_conditions)

evolving = evolving_salinity_mushy_thermodynamics(grid;
    energy_transport,
    salinity_transport = NoSalinityTransport(),
    boundary_conditions)

set!(fixed;
     bulk_salinity = z -> 4 + z,
     temperature = z -> -15 + 3z)

set!(evolving;
     bulk_salinity = z -> 4 + z,
     temperature = z -> -15 + 3z)

initial_salinity = vec(Array(interior(evolving.fields.bulk_salinity)))

for _ in 1:100
    column_energy_time_step!(fixed, 1000)
    column_energy_time_step!(evolving, 1000)
end

diagnostics = (max_internal_energy_difference =
                   maximum(abs.(vec(Array(interior(fixed.fields.internal_energy))) .-
                                vec(Array(interior(evolving.fields.internal_energy))))),
               max_temperature_difference =
                   maximum(abs.(vec(Array(interior(fixed.fields.temperature))) .-
                                vec(Array(interior(evolving.fields.temperature))))),
               max_salinity_drift =
                   maximum(abs.(vec(Array(interior(evolving.fields.bulk_salinity))) .-
                                initial_salinity)))

(diagnostics..., passes = diagnostics.max_internal_energy_difference < 1e-10 &&
                          diagnostics.max_temperature_difference < 1e-10 &&
                          diagnostics.max_salinity_drift == 0)
```

## Stefan Thickness Updates

For a residual interface energy flux ``\delta Q`` that is not retained in the
fixed-grid internal-energy column, the conservative Stefan thickness increment
is

```math
\Delta h = \frac{\Delta t\,\delta Q}{\rho_i L_0}.
```

[`MeltingLimitedSurfaceFlux`](@ref) computes the surface residual by comparing
the requested top flux to the flux required to bring the top cell to
[`complete_melt_energy`](@ref). The applied energy flux warms the column; the
residual is negative for melt and can be passed to
[`column_stefan_thickness_update!`](@ref). The same residual closes
[`column_energy_budget`](@ref) through the `surface_stefan_residual_flux`
keyword.

```@example column_surface_melt
using Oceananigans
using Oceananigans.Fields: set!, interior
using ClimaSeaIce.SeaIceThermodynamics:
    ColumnBoundaryConditions,
    ConductiveTemperatureTransport,
    InsulatingBoundary,
    MeltingLimitedSurfaceFlux,
    QuadraticLiquidusEnergyRelation,
    column_energy_budget,
    column_energy_time_step!,
    column_integrated_energy,
    column_surface_stefan_residual_flux,
    complete_melt_energy,
    internal_energy,
    prescribed_salinity_enthalpy_thermodynamics

grid = RectilinearGrid(size = 1,
                       z = (0, 1),
                       topology = (Flat, Flat, Bounded))

relation = QuadraticLiquidusEnergyRelation(Float64)
S = 5.0
T = -2.0
dt = 100.0
E₀ = internal_energy(relation, T, S)
Ê = complete_melt_energy(relation, S)
excess_flux = 3.0
requested_flux = (Ê - E₀) / dt + excess_flux

thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
    relation,
    salinity_profile = S,
    energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
    boundary_conditions = ColumnBoundaryConditions(
        top = MeltingLimitedSurfaceFlux(flux = requested_flux),
        bottom = InsulatingBoundary()))

set!(thermodynamics; bulk_salinity = S, temperature = T)
initial_energy = column_integrated_energy(thermodynamics)
surface_residual = column_surface_stefan_residual_flux(thermodynamics, dt)
column_energy_time_step!(thermodynamics, dt)
budget = column_energy_budget(thermodynamics, initial_energy, dt;
                              surface_stefan_residual_flux = surface_residual)

abs(first(interior(thermodynamics.fields.internal_energy)) - Ê) < 1e-7 &&
    budget.relative_residual < 1e-12 &&
    surface_residual < 0
```

[`column_stefan_thickness_update!`](@ref) applies this update to a thickness
field, with positive residual flux growing ice and negative residual flux
melting ice. [`column_stefan_thickness_budget`](@ref) reports the corresponding
scalar residual.

```@example column_stefan_update
using ClimaSeaIce.SeaIceThermodynamics:
    QuadraticLiquidusEnergyRelation,
    column_stefan_thickness_budget,
    column_stefan_thickness_change

phase_transitions = QuadraticLiquidusEnergyRelation(Float64).phase_transitions
ρi = 900.0
δQ = 12.0
dt = 3600.0
Δh = column_stefan_thickness_change(phase_transitions, ρi, δQ, dt)
budget = column_stefan_thickness_budget(1.0, 1.0 + Δh,
                                        phase_transitions, ρi, δQ, dt)

budget.relative_residual < 1e-12
```

## Validation

The focused `column_energy` test group checks the thermodynamic relation,
container interfaces, fixed-grid energy budgets, salinity budgets, CPU
allocation discipline, runtime scaling, conservative Stefan thickness updates,
and a manufactured pure-ice conductive mode. The manufactured case verifies
second order spatial convergence of the finite-volume diffusion operator and
first order temporal convergence of the backward-Euler step.
