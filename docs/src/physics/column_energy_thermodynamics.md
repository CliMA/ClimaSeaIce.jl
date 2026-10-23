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

The energy flux closure is

```math
J^E = k \partial_z T + \kappa_E \partial_z E,
```

represented by [`ConductiveTemperatureTransport`](@ref), whose `diffusivity` ``\kappa_E`` defaults to zero. Bulk-salinity transport is either disabled with
[`NoSalinityTransport`](@ref), or stepped with closed boundary scalar diffusion via [`BulkSalinityDiffusion`](@ref).

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

This is the equation assembled at every column step: the right-hand
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

A column is stepped by `SeaIceModel`, which reads the forcing from `model.external_heat_fluxes` and moves the column
interfaces with the Stefan thickness update described below.

```@example column_model
using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: ConductiveTemperatureTransport, IceWaterThermalEquilibrium

grid = RectilinearGrid(size = 8, z = SeaIceColumnDiscretization((0, 1)), topology = (Flat, Flat, Bounded))

thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
    energy_transport = ConductiveTemperatureTransport(conductivity = 2),
    heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = IceWaterThermalEquilibrium(salinity = 0)))

model = SeaIceModel(grid; ice_thermodynamics = thermodynamics, top_heat_flux = 20)

set!(model, h = 1, ℵ = 1)
set!(thermodynamics; temperature = -5)

for _ in 1:24
    time_step!(model, 1hour)
end

first(interior(model.ice_thickness))
```

## Evolving Salinity

The evolving-salinity preset adds `bulk_salinity` to the prognostic fields. With
[`NoSalinityTransport`](@ref), it reduces to the fixed-salinity model for the
same initial salinity field. With [`BulkSalinityDiffusion`](@ref), closed
boundaries conserve column-integrated salinity while reducing salinity variance.

## Fixed and Evolving Presets

The fixed-salinity preset is a strict subcase of the evolving-salinity preset
when salinity transport is disabled. The examples section includes a single-column comparison between a
Bitz-Lipscomb-style fixed-salinity column and an evolving-salinity mushy column
with bulk-salinity diffusion.

## Stefan Thickness Updates

For a residual interface energy flux ``\delta J^E`` that is not retained in the
column internal-energy state, the conservative Stefan thickness increment
is

```math
\Delta h = \frac{\Delta t\,\delta J^E}{\rho_i L_0}.
```

A [`MeltingConstrainedFluxBalance`](@ref) top boundary computes the surface residual by comparing the requested top
flux to the flux required to bring the top cell to [`complete_melt_energy`](@ref). The applied energy flux warms the
column and the residual, negative for melt, moves the surface interface. At a Dirichlet base the residual is the
conductive flux into the base minus the ocean heat flux, and it moves the basal interface.

## Validation

The focused `column_energy` test group checks the thermodynamic relation, container interfaces, stationary and
moving-grid energy budgets, salinity conservation, the melting-limited surface balance, the coupled Stefan update, and a
manufactured pure-ice conductive mode. The manufactured case verifies second order spatial
convergence of the finite-volume diffusion operator and first order temporal
convergence of the backward-Euler step.
