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

On a stationary vertical grid the column solves

```math
\partial_t E = \partial_z J^E + \partial_z I,
\qquad
\partial_t S = \partial_z J_S.
```

The currently implemented energy flux closures are

```math
J^E = k \partial_z T,
\qquad
J^E = \kappa_E \partial_z E,
\qquad
J^E = k \partial_z T + \kappa_E \partial_z E,
```

represented by [`ConductiveTemperatureTransport`](@ref),
[`DiffusiveEnergyTransport`](@ref), and
[`ConductiveAndDiffusiveEnergyTransport`](@ref). Bulk-salinity transport is
either disabled with [`NoSalinityTransport`](@ref), or stepped with closed
boundary scalar diffusion via [`BulkSalinityDiffusion`](@ref).
[`BrineSalinityDiffusion`](@ref) is reserved as an explicit future marker for a
brine-salinity transport equation and is not used by the scalar bulk-salinity
step.

Boundary behavior is configured through the `heat_boundary_conditions = (top, bottom)`
named tuple, while the forcing values are supplied separately as a model-style
`external_heat_fluxes = (top, bottom)` set evaluated through `getflux`.
[`FluxBoundary`](@ref) injects the paired external flux directly across the face,
positive in the increasing vertical-coordinate direction; a resting column with
zero external flux is therefore insulating. [`PrescribedTemperature`](@ref) imposes
a one-sided conductive temperature boundary that is linearized implicitly into the
energy system. At the bottom face this uses

```math
F^E_{1/2} = G_{1/2} (T_1^{n+1} - T_b),
\qquad
G_{1/2} = \frac{2 k_1^n}{\Delta z_1^{n+1}},
```

and the top face uses the corresponding
``F^E_{N+1/2} = G_{N+1/2}(T_t - T_N^{n+1})``. This matches the BL99/Icepack
bottom-ocean-temperature conductance used by the validation replay.
[`MeltingConstrainedFluxBalance`](@ref) imposes a top surface flux but caps the
applied column energy at the complete-melt threshold of the top cell,
returning the excess as a Stefan residual for surface melt.
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
bulk-salinity diffusion step uses the same tridiagonal machinery when bulk
salinity is prognostic. The CPU path uses an allocation-free Thomas sweep over
the same coefficient fields; non-CPU architectures use the Oceananigans batched
tridiagonal solver path.

### Moving Vertical Metric

When the grid uses `MutableVerticalDiscretization`, the vertical coordinate is
treated in conservative moving-coordinate form. Let ``r`` be the reference
coordinate and

```math
z(r, t) = \eta(t) + \mathcal{J}(t) r,
\qquad
\mathcal{J} = \partial z / \partial r,
\qquad
\Delta z_k^n = \mathcal{J}^n \Delta r_k.
```

During one time step the solver treats the conductive/diffusive flux
implicitly and the metric motion explicitly. The implemented thermodynamic
column holds ``\eta`` fixed over the scalar step and represents bottom-fixed
top ablation or top-fixed basal growth by choosing the reference interval and
``\eta`` before updating ``\mathcal{J}``. In the code this metric is Oceananigans'
``\sigma`` field on `MutableVerticalDiscretization`; the previous metric
``\sigma^-`` and current metric ``\sigma^n`` define the face displacement used by
the moving-face flux,

```math
\delta z_{g,k+1/2}
= z_{k+1/2}^{n+1} - z_{k+1/2}^n
= \left(\mathcal{J}^{n+1} - \mathcal{J}^n\right) r_{k+1/2}
= \left(\sigma^n_{k+1/2} - \sigma^-_{k+1/2}\right) r_{k+1/2}.
```

With the stationary-grid conductive/diffusive flux renamed ``F^E`` to avoid
confusing it with the Jacobian ``\mathcal{J}``, the moving-coordinate equation
represented by the discretization is

```math
\partial_t(\mathcal{J} E)
= \partial_r(F^E + I) + \partial_r(\dot z_g E^{up}),
\qquad
\partial_t(\mathcal{J} S)
= \partial_r F^S + \partial_r(\dot z_g S^{up}),
```

where ``\dot z_g`` is the grid-face velocity and ``E^{up}`` and ``S^{up}``
are the upwind cell-centered values swept by each moving face. The
finite-volume energy update used by the solver is therefore

```math
\mathcal{J}^{n+1} \Delta r_k E_k^{n+1}
= \mathcal{J}^n \Delta r_k E_k^n
+ \left[
\delta z_{g,k+1/2} E^{up}_{k+1/2}
- \delta z_{g,k-1/2} E^{up}_{k-1/2}
\right]
+ \Delta t \left[
(F^E + I)_{k+1/2}^{n+1}
- (F^E + I)_{k-1/2}^{n+1}
\right],
```

where ``E^{up}`` is the piecewise-constant cell value swept across the moving
face, ``F^E`` is the conductive/diffusive internal-energy flux (the same
stationary-grid flux denoted ``J^E`` above), and ``I`` is the shortwave flux.
Thus, after division by the current physical layer thickness,

```math
E_k^{n+1}
= \frac{\mathcal{J}^n}{\mathcal{J}^{n+1}} E_k^n
+ \frac{\Delta t}{\Delta z_k^{n+1}}
\left[(F^E + I)_{k+1/2}^{n+1}
- (F^E + I)_{k-1/2}^{n+1}\right]
+ \frac{
\delta z_{g,k+1/2} E^{up}_{k+1/2}
- \delta z_{g,k-1/2} E^{up}_{k-1/2}}
{\Delta z_k^{n+1}}.
```

This is the equation assembled by `column_energy_time_step!`: the right-hand
side contains the old concentration scaled by the metric ratio
``\sigma^-_k / \sigma^n_k``, the explicit swept-face enthalpy integral divided
by ``\Delta z_k^{n+1}``, and the implicit conductive/diffusive flux divergence
using current physical distances. The tridiagonal solve therefore advances
``E^{n+1}`` as a concentration while conserving the layer integral
``\mathcal{J}^{n+1}\Delta r_k E_k^{n+1}``.

The implicit diffusion coefficients use current physical distances,

```math
a_{k,k+1}
= \frac{\Delta t\,D_{k+1/2}}
       {\Delta z_k^{n+1}\,\Delta z_{k+1/2}^{n+1}},
```

where ``D`` is the effective energy diffusivity from the linearized temperature
transport. Prognostic bulk salinity uses the analogous conservative equation

```math
\mathcal{J}^{n+1} \Delta r_k S_k^{n+1}
= \mathcal{J}^n \Delta r_k S_k^n
+ \left[
\delta z_{g,k+1/2} S^{up}_{k+1/2}
- \delta z_{g,k-1/2} S^{up}_{k-1/2}
\right]
+ \Delta t \left[F^S_{k+1/2}^{n+1} - F^S_{k-1/2}^{n+1}\right].
```

Thus a no-flux moving-boundary step preserves uniform energy and salinity
concentrations while changing the layer integrals with the column thickness. A
nonuniform no-flux step uses the cell swept by the moving face: for
``\Delta z_g > 0`` the face samples the cell above, and for
``\Delta z_g < 0`` it samples the cell below. This is equivalent to a
piecewise-constant conservative overlap remap when a face crosses at most one
cell during the step. Boundary faces default to the adjacent interior value,
which preserves uniform concentrations during no-flux expansion. When growth
creates material with a distinct enthalpy, the volumetric internal energy swept
in by the moving boundary can be prescribed separately while retaining the same
imposed boundary flux. When
``\mathcal{J}^{n+1}=\mathcal{J}^n``, the moving-face term vanishes and these
equations reduce exactly to the stationary-grid system above.

```@example column_energy_step
using Oceananigans
using Oceananigans.Fields: set!
using ClimaSeaIce.SeaIceThermodynamics:
    ConductiveTemperatureTransport,
    FluxBoundary,
    prescribed_salinity_enthalpy_thermodynamics,
    column_energy_budget,
    column_energy_time_step!,
    column_integrated_energy
using ClimaSeaIce: SeaIceColumnDiscretization

# Forcing lives in a model-style `external_heat_fluxes = (top, bottom)` set evaluated through `getflux`,
# decoupled from the `heat_boundary_conditions` behavior. These scalar fluxes ignore the clock and the
# (empty) model fields.
clock = Clock(time = 0.0)
model_fields = NamedTuple()
heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary())
external_heat_fluxes = (top = 5.0, bottom = 2.0)

grid = RectilinearGrid(size = 8,
                       z = SeaIceColumnDiscretization((0, 1)),
                       topology = (Flat, Flat, Bounded))

thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
    salinity_profile = 0.0,
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
    heat_boundary_conditions)

set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> -12 + 4z)
initial_energy = column_integrated_energy(thermodynamics)
column_energy_time_step!(thermodynamics, external_heat_fluxes, clock, model_fields, 5e3)
budget = column_energy_budget(thermodynamics, external_heat_fluxes, clock, model_fields, initial_energy, 5e3)

budget.relative_residual < 1e-11
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
    heat_boundary_conditions)

shortwave_fluxes = (top = 0.0, bottom = 0.0)
set!(shortwave_thermodynamics; bulk_salinity = 0.0, temperature = -10.0)
initial_energy = column_integrated_energy(shortwave_thermodynamics)
column_energy_time_step!(shortwave_thermodynamics, shortwave_fluxes, clock, model_fields, 100.0)
budget = column_energy_budget(shortwave_thermodynamics, shortwave_fluxes, clock, model_fields, initial_energy, 100.0)

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
    ConductiveTemperatureTransport,
    FluxBoundary,
    evolving_salinity_mushy_thermodynamics,
    column_integrated_salinity,
    column_salt_budget,
    column_salinity_time_step!
using ClimaSeaIce: SeaIceColumnDiscretization

grid = RectilinearGrid(size = 8,
                       z = SeaIceColumnDiscretization((0, 1)),
                       topology = (Flat, Flat, Bounded))

thermodynamics = evolving_salinity_mushy_thermodynamics(grid;
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
    salinity_transport = BulkSalinityDiffusion(diffusivity = 1e-4),
    heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary()))

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
The examples section also includes a single-column comparison between a
Bitz-Lipscomb-style fixed-salinity column and an evolving-salinity mushy column
with bulk-salinity diffusion.

```@example column_equivalence
using Oceananigans
using Oceananigans.Fields: set!, interior
using ClimaSeaIce.SeaIceThermodynamics:
    ConductiveTemperatureTransport,
    FluxBoundary,
    NoSalinityTransport,
    column_energy_time_step!,
    evolving_salinity_mushy_thermodynamics,
    prescribed_salinity_enthalpy_thermodynamics
using ClimaSeaIce: SeaIceColumnDiscretization

clock = Clock(time = 0.0)
model_fields = NamedTuple()

grid = RectilinearGrid(size = 8,
                       z = SeaIceColumnDiscretization((0, 1)),
                       topology = (Flat, Flat, Bounded))

heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary())
external_heat_fluxes = (top = 0.0, bottom = 0.0)
energy_transport = ConductiveTemperatureTransport(conductivity = 2.0)

fixed = prescribed_salinity_enthalpy_thermodynamics(grid;
    salinity_profile = 0.0,
    energy_transport,
    heat_boundary_conditions)

evolving = evolving_salinity_mushy_thermodynamics(grid;
    energy_transport,
    salinity_transport = NoSalinityTransport(),
    heat_boundary_conditions)

set!(fixed;
     bulk_salinity = z -> 4 + z,
     temperature = z -> -15 + 3z)

set!(evolving;
     bulk_salinity = z -> 4 + z,
     temperature = z -> -15 + 3z)

initial_salinity = vec(Array(interior(evolving.fields.bulk_salinity)))

for _ in 1:100
    column_energy_time_step!(fixed, external_heat_fluxes, clock, model_fields, 1000)
    column_energy_time_step!(evolving, external_heat_fluxes, clock, model_fields, 1000)
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

For a residual interface energy flux ``\delta J^E`` that is not retained in the
column internal-energy state, the conservative Stefan thickness increment
is

```math
\Delta h = \frac{\Delta t\,\delta J^E}{\rho_i L_0}.
```

A [`MeltingConstrainedFluxBalance`](@ref) top boundary computes the surface
residual by comparing the requested top flux to the flux required to bring the
top cell to [`complete_melt_energy`](@ref). The applied energy flux warms the
column; the residual is negative for melt and can be passed to
[`column_stefan_thickness_update!`](@ref). The same residual closes
[`column_energy_budget`](@ref) through the `surface_stefan_residual_flux`
keyword, and is computed by
[`compute_column_surface_stefan_residual_flux!`](@ref).

```@example column_surface_melt
using Oceananigans
using Oceananigans.Fields: set!, interior
using ClimaSeaIce.SeaIceThermodynamics:
    ConductiveTemperatureTransport,
    FluxBoundary,
    MeltingConstrainedFluxBalance,
    QuadraticLiquidusEnergyRelation,
    column_energy_budget,
    column_energy_time_step!,
    column_integrated_energy,
    compute_column_surface_stefan_residual_flux!,
    complete_melt_energy,
    internal_energy,
    prescribed_salinity_enthalpy_thermodynamics
using ClimaSeaIce: SeaIceColumnDiscretization

clock = Clock(time = 0.0)
model_fields = NamedTuple()

grid = RectilinearGrid(size = 1,
                       z = SeaIceColumnDiscretization((0, 1)),
                       topology = (Flat, Flat, Bounded))

relation = QuadraticLiquidusEnergyRelation(Float64)
S = 5.0
T = -2.0
dt = 100.0
E₀ = internal_energy(relation, T, S)
Ê = complete_melt_energy(relation, S)
excess_flux = 3.0
requested_flux = (Ê - E₀) / dt + excess_flux

# Upward-positive convention: a warming surface flux (heat into the ice) is negative.
external_heat_fluxes = (top = -requested_flux, bottom = 0.0)

thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
    relation,
    salinity_profile = S,
    energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
    heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(),
                                bottom = FluxBoundary()))

set!(thermodynamics; bulk_salinity = S, temperature = T)
initial_energy = column_integrated_energy(thermodynamics)

residual_flux = Field{Center, Center, Nothing}(grid)
compute_column_surface_stefan_residual_flux!(residual_flux, thermodynamics, external_heat_fluxes, clock, model_fields, dt)
surface_residual = first(interior(residual_flux))

column_energy_time_step!(thermodynamics, external_heat_fluxes, clock, model_fields, dt)
budget = column_energy_budget(thermodynamics, external_heat_fluxes, clock, model_fields, initial_energy, dt;
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

## Split Thickness Remapping

CICE/Icepack `thickness_changes` first advances the temperature solve, then
repartitions layer enthalpy onto the new equal-layer geometry. ClimaSeaIce
exposes the corresponding split operation with
[`column_energy_thickness_remap!`](@ref). The helper conservatively remaps
layer-averaged internal energy from `source_faces` to `target_faces`, fills
uncovered target intervals with a supplied `fill_energy`, and recomputes
diagnostics. Passing `bulk_salinity` resets fixed-salinity profiles after the
remap instead of remapping salinity as a prognostic tracer.

For old source intervals ``[z^-_{\ell-1/2}, z^-_{\ell+1/2}]`` and new target
intervals ``[z^+_{k-1/2}, z^+_{k+1/2}]``, the remapped internal energy is

```math
E^+_k =
\frac{1}{\Delta z^+_k}
\left[
\sum_\ell E^-_\ell
\left|[z^-_{\ell-1/2}, z^-_{\ell+1/2}]
      \cap
      [z^+_{k-1/2}, z^+_{k+1/2}]\right|
+ E_\mathrm{fill}
\left(\Delta z^+_k -
\sum_\ell
\left|[z^-_{\ell-1/2}, z^-_{\ell+1/2}]
      \cap
      [z^+_{k-1/2}, z^+_{k+1/2}]\right|
\right)
\right].
```

This is the finite-volume counterpart of the moving-metric swept-face term
above, applied as a split geometry update. It is used when the CICE-compatible
sequence requires a temperature solve on the old equal-layer grid followed by
Icepack-style equal-layer repartitioning on the new thickness.

For top ablation, `source_faces` and `target_faces` share the bottom face and
the top part of the old column is removed. For basal growth, the old source
faces are offset upward by the growth amount and the exposed lower interval is
filled with the new-ice enthalpy.

```@example column_thickness_remap
using Oceananigans
using Oceananigans.Fields: set!, interior
using ClimaSeaIce.SeaIceThermodynamics:
    ConductiveTemperatureTransport,
    FluxBoundary,
    QuadraticLiquidusEnergyRelation,
    column_energy_thickness_remap!,
    conservative_column_remap,
    prescribed_salinity_enthalpy_thermodynamics
using ClimaSeaIce: SeaIceColumnDiscretization

grid = RectilinearGrid(size = 4,
                       z = SeaIceColumnDiscretization((0, 1)),
                       topology = (Flat, Flat, Bounded))

relation = QuadraticLiquidusEnergyRelation(Float64)
thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
    relation,
    salinity_profile = 0,
    energy_transport = ConductiveTemperatureTransport(conductivity = 0),
    heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary()))

set!(thermodynamics; bulk_salinity = 0, temperature = z -> -12 + 8z)

source_faces = collect(range(0, 4; length = 5))
target_faces = collect(range(0, 3; length = 5))
source_energy = vec(Array(interior(thermodynamics.fields.internal_energy)))
expected_energy = conservative_column_remap(source_energy,
                                            source_faces,
                                            target_faces)

column_energy_thickness_remap!(thermodynamics,
                               source_faces,
                               target_faces;
                               bulk_salinity = z -> 1 + z)

maximum(abs.(vec(Array(interior(thermodynamics.fields.internal_energy))) .-
             expected_energy)) < 1e-12
```

## Validation

The focused `column_energy` test group checks the thermodynamic relation,
container interfaces, stationary and moving-grid energy budgets, salinity
budgets, CPU allocation discipline, runtime scaling, conservative Stefan
thickness updates, conservative split thickness remapping, and a manufactured
pure-ice conductive mode. The manufactured case verifies second order spatial
convergence of the finite-volume diffusion operator and first order temporal
convergence of the backward-Euler step.
