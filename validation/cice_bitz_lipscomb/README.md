# CICE Bitz-Lipscomb Column Validation

This directory records the CICE/Icepack reference setup for validating the
ClimaSeaIce fixed-salinity column-energy thermodynamics against CICE column
mode with Icepack BL99 thermodynamics.

Status: pre-CICE proof gates are in place and local no-snow CICE column cases
now build, run, parse instantaneous CICE history, map parsed states into the
ClimaSeaIce fixed-salinity column relation, and record reset-each-interval
fixed-grid forced replays. The column energy solve now has a
`MutableVerticalDiscretization` path that evolves energy and prognostic salinity
with the moving vertical Jacobian and explicit moving-face enthalpy/salinity
flux. The CICE comparison scripts now report both the older fixed-grid/remap
diagnostics and a moving-metric replay. The strict BL99 validation gate now
passes for all eight no-snow histories by using the source-level Icepack BL99
temperature matrix, CICE-style thickness-change remapping, and the source-level
prognostic-thickness replay. The current
`cold_conductive_relaxation` CICE runs use a validation-only
balanced-bottom-flux hook and pass the fixed-thickness
transition gate for both `MU71` and `bubbly`; the fixed-thickness replay now
holds the predicted thickness fixed while still reporting CICE's observed `hi`
drift separately. Matching the Icepack BL99
`cp_ice = 2106` and `cp_ocn = 4218` constants brings fixed-cold, basal-growth,
and controlled surface-warming reset temperature, per-level shifted internal
energy, and column-energy diagnostics below 1%. A ClimaSeaIce source-level
Icepack BL99 temperature-matrix replay followed by CICE-style conservative
thickness remapping now passes all thickness-changing reset diagnostics,
including forced surface ablation. The forced surface-ablation cases produce
the intended top melt for both conductivity variants; their raw fixed-grid
reset replays miss the 1% target, and ClimaSeaIce fixed-grid remapping still
leaves per-level shifted internal-energy diagnostics above 1% for ablation.
The general ClimaSeaIce CICE-thickness-forced replay passes fixed-cold,
controlled warming, and basal growth by the official temperature gate, but
forced ablation still fails in the semi-implicit split-solve path. A sequential
source-level Icepack matrix replay with the CICE thickness history passes the
official temperature, CICE-negative-enthalpy, and column-energy gates for all
histories. The moving-metric replay confirms that the conservative
`MutableVerticalDiscretization` Jacobian equations are necessary, but CICE
agreement for growth and ablation still comes from Icepack's split
`thickness_changes` enthalpy redistribution rather than from a single
continuous moving-grid diffusion solve.

## Reference Checkout

Local source checkout:

```text
/private/tmp/cice-validation/CICE
```

Reference commits inspected:

| component | commit |
| --- | --- |
| CICE | `3c2064659ccf693e644e0cfd21140565ddb8e355` |
| Icepack submodule | `68097c59aceb2040f8c0cef189bccfe7590a8a32` |

## Source References

The source line numbers below refer to the checkout above.

| Topic | Source reference | Notes for ClimaSeaIce matching |
| --- | --- | --- |
| BL99 constants | `icepack/columnphysics/icepack_parameters.F90:68`, `:94`, `:106-109`, `:120-125`, `:142`, `:148` | `puny = 1e-11`, `rhoi = 917`, `cp_ice = 2106`, `cp_ocn = 4218`, `depressT = 0.054`, `kice = 2.03`, `saltmax = 3.2`, default `ktherm = 1`, default `conduct = 'bubbly'`. |
| Salinity profile setup | `icepack/columnphysics/icepack_therm_shared.F90:233-268` | `icepack_init_salinity` fills layer-center salinity and sets `l_brine` from `saltmax`. |
| Salinity profile formula | `icepack/columnphysics/icepack_therm_shared.F90:285-302` | `S(zn) = saltmax / 2 * (1 - cos(pi * zn^(0.407 / (0.573 + zn))))`; `zn` is normalized downward depth from the top surface. |
| Initial BL99 enthalpy | `icepack/columnphysics/icepack_therm_shared.F90:344-355` | For `ktherm != 2`, CICE stores negative enthalpy `qin = -rhoi * (cp_ice * (Tmlt - Ti) + Lfresh * (1 - Tmlt / Ti) - cp_ocn * Tmlt)`. |
| Temperature inversion | `icepack/columnphysics/icepack_therm_shared.F90:68-99` | `calculate_Tin_from_qin` solves a quadratic and clamps `Tin <= Tmltk`. |
| Runtime BL99 enthalpy update | `icepack/columnphysics/icepack_therm_bl99.F90:630-638` | Recomputes `zqin` from `zTin`, `Tmlts`, and the BL99 brine-pocket expression. |
| BL99 heat capacity iteration | `icepack/columnphysics/icepack_therm_bl99.F90:343-357` | Iteration uses `ci = cp_ice - Lfresh * Tmlt / (zTin * Tin_init)` when brine pockets are active. |
| Liquidus temperature | `icepack/columnphysics/icepack_therm_bl99.F90:236-240` | Runtime layer melting temperature is `Tmlts(k) = -zSin(k) * depressT`. |
| Conductivity laws | `icepack/columnphysics/icepack_therm_bl99.F90:824-875` | `conduct = 'MU71'` and `conduct = 'bubbly'` both use `min(-puny, zTin)` in the salinity denominator and clamp with `kimin = 0.10`. |
| Interface conductance | `icepack/columnphysics/icepack_therm_bl99.F90:877-908` | Interior ice `kh = 2 * K_{k-1} * K_k / ((K_{k-1} + K_k) * hilyr)`, in conductance units `W m^-2 K^-1`. ClimaSeaIce face conductivity must multiply this by the center-to-center spacing for its grid-operator form. |
| Surface temperature constraint | `icepack/columnphysics/icepack_therm_bl99.F90:387-399`, `:459-479`, `:509-548`, `:657-674` | Iteration enforces `Tsf <= 0 C`; when `Tsf = 0 C`, convergence also requires surface flux not to exceed conductive capacity. |
| Interior melt limiting | `icepack/columnphysics/icepack_therm_bl99.F90:591-599`, `:697-703` | Layers exceeding melting temperature are reset to `Tmlt`, excess energy is tracked, and interface conductance can be reduced next iteration. |
| Bottom conductive flux | `icepack/columnphysics/icepack_therm_bl99.F90:683-687` | Bottom flux uses the bottom interface conductance and includes any extra energy routed out of ice. |
| Energy residual check | `icepack/columnphysics/icepack_therm_bl99.F90:678-690` | Iteration checks new minus initial internal energy against top conduction, bottom conduction, and absorbed shortwave. |
| CICE state initialization | `cicecore/cicedyn/general/ice_init.F90:3639-3650` | CICE initializes `qin` and stores `nt_qice` and `nt_sice` layer tracers from the Icepack profile. |
| Enthalpy helper used by BGC remap | `icepack/columnphysics/icepack_zbgc_shared.F90:216-230` | Same BL99 negative enthalpy expression appears in `calculate_qin_from_Sin`. |
| Equal-layer remap helper | `icepack/columnphysics/icepack_therm_shared.F90:485-568` | `adjust_enthalpy` computes overlap integrals from old unequal layers to new equal-thickness layers for enthalpy-like layer quantities. |
| Top-ablation thickness update | `icepack/columnphysics/icepack_therm_vertical.F90:1488-1515`, `:1718-1738` | For `ktherm != 2`, top melt removes ice using `qm = zqin`, then calls `adjust_enthalpy` to repartition remaining enthalpy onto equal layers. |
| Basal congelation thickness update | `icepack/columnphysics/icepack_therm_vertical.F90:1345-1401`, `:1718-1738` | For `ktherm != 2`, `thickness_changes` computes `qbot` from bottom salinity and `Tbot`, grows the bottom layer, then calls `adjust_enthalpy` to repartition to equal layers. |
| Vertical tracer remap for added basal ice | `icepack/columnphysics/icepack_therm_itd.F90:800-868` | `update_vertical_tracers` handles bottom-added ice by overlap remapping in the mushy `ktherm = 2` vertical tracer path. |
| BL99 open-water new-ice enthalpy/salinity path | `icepack/columnphysics/icepack_therm_itd.F90:1452-1475`, `:1606-1688` | For `ktherm != 2`, open-water/frazil new ice uses `qi0new = -rhoi * Lfresh`, the BL99 salinity profile, and distribution through the existing column; this is separate from the basal congelation path active in the winter growth cases here. |
| Calm ocean forcing used by the validation templates | `cicecore/cicedyn/general/ice_flux.F90:1047-1060`, `cicecore/cicedyn/general/ice_forcing.F90:4162-4178` | The column setup initializes `frzmlt = 0`, `sst = Tf`, and `qdp = 0`. With `oceanmixed_ice = .false.`, there is no mixed-layer heat source to balance a cold bottom conductive flux. |
| Standalone mixed-layer gate | `cicecore/drivers/standalone/cice/CICE_RunMod.F90:490-497` | The standalone driver calls `ocean_mixed_layer` only when `oceanmixed_ice` is enabled; the validation templates currently disable it. |
| Fixed-thickness cold validation hook | `validation/cice_bitz_lipscomb/run_cice_cases.sh` | With `CICE_BALANCE_COLD_BOTTOM_FLUX=1`, the runner injects a small Icepack source hook after the BL99 temperature solve. It feeds the computed `fcondbotn` back to `thickness_changes` as the bottom ocean heat flux, preserving the temperature solve while preventing basal growth. |

## Local Proof Gates

The ClimaSeaIce proof gates currently covered by `TEST_GROUP=column_energy`
include:

- BL99 salinity midpoint values, endpoints, and coordinate reversal.
- Fixed-salinity brine-pocket temperature inversion, pure-ice limit,
  complete-melt threshold, finite-difference derivatives, and heat-capacity
  identity.
- `MaykutUntersteinerConductivity` and `BubblyBrineConductivity` scalar
  formulas, limiter behavior, and harmonic face conductivities.
- `MutableVerticalDiscretization` metric integration, including moving-face
  enthalpy/salinity fluxes, for no-flux energy, boundary-flux energy, and
  no-transport prognostic salinity steps. The energy proof gate now also covers
  moving-boundary fill enthalpy, which is needed for basal growth where new ice
  enters with the Icepack bottom-growth enthalpy rather than the old bottom-cell
  value.
- Conservative piecewise-constant layer remapping for layer-integrated column
  quantities, including top-ablation and basal-growth fill-value tests. The
  `ColumnEnergyThermodynamics` remap helper now resets fixed salinity from a
  supplied normalized-depth profile after remapping, matching the BL99
  validation requirement that salinity is recomputed rather than transported.

The public ClimaSeaIce names used for this BL99 validation are
`FixedDrainedIceSalinityProfile`, `FixedSalinityBrinePocketEnergyRelation`,
`MaykutUntersteinerConductivity`, `BubblyBrineConductivity`,
`MeltingConstrainedSurfaceFluxBalance`, and
`OceanFreezingTemperatureBoundary`. The last two are descriptive aliases for
the existing `MeltingConstrainedFluxBalance` and `IceWaterThermalEquilibrium`
boundary conditions.

The latest lightweight gate status is in
`validation/cice_bitz_lipscomb/results/pre_cice_dashboard.csv`.

The one-figure visual summary is generated by
`validation/cice_bitz_lipscomb/plot_cice_single_column_validation.jl` and
written to
`validation/cice_bitz_lipscomb/results/cice_single_column_validation.png`.

The current multi-case CICE evidence is aggregated in
`validation/cice_bitz_lipscomb/results/cice_case_comparison_summary.csv`:

| case | conduct | transition | reset replay | Icepack-matrix remapped reset | thickness-remapped reset | CICE-thickness-forced replay | free replay | note |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `cold_conductive_relaxation` | `MU71` | pass | pass, `0.00103279653852` | pass, `5.14638749576e-06` | pass, `0.00103304883639` | pass, `0.00799516820487` | pass, `0.00799365465354` | Balanced-bottom fixed-thickness run; `delta_hi = 3.0e-6 m`, `delta(hi/aice) = 3.30000099003e-6 m`, no basal growth. |
| `cold_conductive_relaxation` | `bubbly` | pass | pass, `0.00110474853209` | pass, `5.74253411496e-07` | pass, `0.00110506393058` | pass, `0.00778523135678` | pass, `0.0077841825518` | Balanced-bottom fixed-thickness run; `delta_hi = -1.2e-6 m`, `delta(hi/aice) = -4.00000319978e-7 m`, no basal growth. |
| `basal_growth` | `MU71` | pass | pass, `0.00263944227201` | pass, `9.26538101245e-06` | pass, `0.00107403724222` | pass, `0.00875578804266` | fail, `0.0485297179925` | Uses the unbalanced winter bottom-flux setup to retain measurable basal congelation. |
| `basal_growth` | `bubbly` | pass | pass, `0.00283466651198` | pass, `1.06792979317e-05` | pass, `0.00113305358948` | pass, `0.00852913566164` | fail, `0.050180069487` | Same basal-growth setup with CICE/Icepack `conduct = 'bubbly'`. |
| `surface_warming` | `MU71` | pass | pass, `0.0012871994515` | pass, `3.84029189641e-06` | pass, `0.000357230763411` | pass, `0.0022686519764` | fail, `0.0205210326828` | Prescribed `8.0 W m^-2` surface flux warms the column while limiting ice loss to `0.0045599 m` over two days. |
| `surface_warming` | `bubbly` | pass | pass, `0.0012921956637` | pass, `3.89723767928e-06` | pass, `0.000363279184033` | pass, `0.00221754423547` | fail, `0.0200419350573` | Same controlled warming setup with CICE/Icepack `conduct = 'bubbly'`. |
| `surface_ablation` | `MU71` | pass | fail, `0.0190483763055` | pass, `0.000919646651979` | pass, `0.0077718269457` | fail, `0.0542395682928` | fail, `0.43770623144` | Prescribed positive surface flux removes `0.10343 m` grid-cell ice volume; `hi/aice` decreases by `0.0531144458597 m`. |
| `surface_ablation` | `bubbly` | pass | fail, `0.0191031535505` | pass, `0.000917650982306` | pass, `0.007828341593` | fail, `0.0523392114103` | fail, `0.421717149327` | Same forced ablation setup; `hi/aice` decreases by `0.0531160845177 m`. |

The aggregate also includes thickness-rate audit columns:
`integrated_congel_m`, `integrated_meltt_m`, `integrated_meltb_m`,
`rate_integrated_delta_hi_m`, `thickness_rate_budget_residual_m`,
`thickness_rate_budget_status`, `initial_aice`, `final_aice`, `delta_aice`,
`initial_thermodynamic_hi_m`, `final_thermodynamic_hi_m`,
`delta_thermodynamic_hi_m`, `prognostic_thickness_replay_status`, and
`strict_bl99_validation_status`. The `hi` columns are CICE grid-cell ice volume
per horizontal area, while `thermodynamic_hi_m = hi / aice` diagnoses the
single-category thermodynamic thickness associated with `Tinz`. The CICE
rate-integrated thickness budgets and `strict_bl99_validation_status` pass for
all eight histories; the largest residual is about `3.24e-4 m` in the
forced surface-ablation cases. The integrated top melt is `0.103046700238 m`
for `surface_ablation` with `MU71` and `0.103035348451 m` with `bubbly`; the
integrated basal growth is `0.0518581154167 m` for `basal_growth` with `MU71`
and `0.056307245 m` with `bubbly`.

The aggregate now also reports a first prognostic-thickness forced replay. This
diagnostic applies the CICE-compatible top flux split, routes the
de-area-weighted top-flux residual into surface melt, computes basal growth
from the predicted bottom conductive residual, converts residual energy to
thickness with the BL99 negative-enthalpy magnitude, and remaps energy onto the
predicted thickness. The CICE history fields `fsurf_ai` and `fcondtop_ai` are
weighted by ice area, so the analyzer divides them by parsed `aice` before
using them as per-ice-area column fluxes. It passes the prognostic
thickness/thermal gate for controlled warming (`0.228143603683` for `MU71`,
`0.223090386649` for `bubbly`), basal growth (`0.88413961273`,
`0.86122201312`), and the fixed-cold no-growth replay (`0.799362278917`,
`0.778412457441`) once fixed-thickness intervals are held fixed in the
predicted replay. It still fails the forced-ablation temperature gate
(`3.57820733413`, `3.49785711177`). For forced ablation the predicted top melt is
`0.102985686902 m` for `MU71` and `0.102975950388 m` for `bubbly`, within the
1%/1 mm thickness gate against CICE. The CICE-observed fixed-cold `delta_hi`
values remain reported separately (`3.0e-6 m` for `MU71` and `-1.2e-6 m` for
`bubbly`), but they are not applied as Stefan growth/melt in this
fixed-thickness replay.

The same aggregate reports `icepack_matrix_prognostic_thickness_forced_*`, a
source-level BL99 temperature-matrix version of the prognostic-thickness replay.
It passes controlled warming (`0.0292835280463` for `MU71`,
`0.0286184463132` for `bubbly`) and basal growth (`0.00504702760601`,
`0.00519788344341`), fixed-cold (`0.000370455930163` and
`0.000196801198451`), and, after de-area-weighting the CICE top fluxes, forced
ablation (`0.990152772454` and `0.983805215958`). The source-level
fixed-thickness drift is zero in the predicted replay. The source-level
forced-ablation top melt is `0.10303560405 m` for `MU71` and
`0.103024427872 m` for `bubbly`, only about `1.1e-5 m` below the integrated
CICE `meltt` history in each conductivity variant. This source-level
prognostic replay is the official strict prognostic-thickness gate used by
`strict_bl99_validation_status`.

The aggregate CSV also reports shifted BL99 internal-energy, CICE
negative-enthalpy, and column-energy errors. The shifted internal-energy
relative floor is frozen from the minimum nonzero reconstructed `|E|` in the
cold conductive MU71 history: `2.37722677947e6 J m^-3`. The CICE
negative-enthalpy floor is frozen from the same history at
`2.95867899346e8 J m^-3`, using `q_i = E - rho_i L0`. With the shifted floor,
representative remapped reset diagnostics are:

| case | conduct | remapped temperature | remapped internal energy | remapped column energy |
| --- | --- | --- | --- | --- |
| `cold_conductive_relaxation` | `MU71` | pass, `0.00103304883639` | pass, `0.00871786496474` | pass, `1.29116797151e-05` |
| `cold_conductive_relaxation` | `bubbly` | pass, `0.00110506393058` | pass, `0.00929842173781` | pass, `1.45141432524e-05` |
| `basal_growth` | `MU71` | pass, `0.00107403724222` | pass, `0.00912957581408` | pass, `1.4564561509e-05` |
| `basal_growth` | `bubbly` | pass, `0.00113305358948` | pass, `0.00963171963647` | pass, `1.64095143675e-05` |
| `surface_warming` | `MU71` | pass, `0.000357230763411` | pass, `0.00089974495413` | pass, `3.13040274025e-07` |
| `surface_warming` | `bubbly` | pass, `0.000363279184033` | pass, `0.000915045669336` | pass, `3.18420119917e-07` |
| `surface_ablation` | `MU71` | pass, `0.0077718269457` | fail, `0.0195520787001` | pass, `1.94602446253e-06` |
| `surface_ablation` | `bubbly` | pass, `0.007828341593` | fail, `0.0196562450921` | pass, `1.97892721208e-06` |

CICE negative enthalpy passes for the same remapped ablation diagnostics:
`0.000152277126163` for `MU71` and `0.000153054745318` for `bubbly`.

Temperature acceptance ratios use the validation gate directly: relative
temperature error divided by `0.01`, except layers with `|T_CICE| < 2 C` may
instead pass by absolute error divided by `0.02 C`. The aggregate reports
paired `*_temperature_acceptance_status` and
`*_temperature_acceptance_ratio` columns for the reset, remapped reset,
CICE-thickness-forced, moving-metric, source-level matrix, top-ablation reset,
and free-running replay paths.

The aggregate also includes `icepack_temperature_matrix_reset_*` columns. This
source-level diagnostic uses the Icepack BL99 temperature matrix with CICE's
top conductive flux and bottom ocean temperature from each history interval. It
reproduces the fixed-cold CICE temperature profiles to `7.96658313105e-07` for
`MU71` and `1.01110189246e-06` for `bubbly`; shifted internal-energy errors are
`1.1539840022e-06` and `2.32760063197e-06`.

The aggregate additionally reports
`icepack_temperature_matrix_thickness_remapped_reset_*`, which applies the
source-level temperature solve and then the CICE-style conservative enthalpy
remap used by `thickness_changes`. This source-level remapped diagnostic passes
all eight histories. For the forced ablation cases, maximum relative
temperature errors are `0.000919646651979` for `MU71` and
`0.000917650982306` for `bubbly`; shifted internal-energy errors are
`0.002322920769` and `0.0023133988512`.

The aggregate also reports
`icepack_temperature_matrix_cice_thickness_forced_*`, a sequential diagnostic
that uses the source-level BL99 temperature matrix while following the CICE
thickness history. This passes the official temperature gate for all eight
histories. The relative-only 1% status still fails forced ablation: maximum
relative temperature errors are `0.0137447491872` for `MU71` and
`0.0134220179813` for `bubbly`, but the corresponding temperature acceptance
ratios are `0.877625726166` and `0.866657092284`. Shifted per-level
internal-energy errors still fail ablation at `0.0344588185143` and
`0.0335931548341`; the same profiles pass the CICE negative-enthalpy metric at
`0.000285545531375` and `0.000293443719422`, and column-energy errors remain
small at `7.77190706404e-05` and `6.08825826508e-05`.

The general ClimaSeaIce CICE-thickness-forced path now uses a prescribed bottom
ocean temperature in the semi-implicit energy system. It passes fixed-cold,
controlled-warming, and basal-growth temperature gates (`0.875578804266` for
`MU71` growth and `0.852913566164` for `bubbly` growth), while forced ablation
still fails (`3.46328913403` and `3.37953270785`) even though its CICE
negative-enthalpy errors are below 1%. This separates the remaining ablation
mismatch from the bottom boundary, source-level BL99 matrix, and CICE enthalpy
reconstruction.

The aggregate also reports `moving_metric_cice_thickness_forced_*`, which runs
the general ClimaSeaIce step on `MutableVerticalDiscretization` with the CICE
`hi` sequence represented by the moving Jacobian and upwind moving-face
enthalpy flux. It passes the fixed-cold histories by temperature
(`0.00799353060986` for `MU71` and `0.00778523000377` for `bubbly`) and passes
the controlled-warming histories by temperature and shifted internal energy
(`0.00245073224801` and `0.00240268783156` temperature error). Basal growth
still fails the official gate (`7.66561191628` and `7.98031241752`) because
the current mutable-grid diagnostic does not yet use the CICE bottom-growth
fill enthalpy and equal-layer repartitioning. Forced ablation is closer after
correcting the swept-face upwind direction, but still fails the official gate
(`3.72528155254` and `3.6419695233`), confirming that top-ablation agreement
needs the full Icepack `thickness_changes` enthalpy redistribution in addition
to the moving metric.
A top-fixed continuous moving-grid growth experiment with the correct basal
fill enthalpy performed worse than the reported bottom-fixed metric diagnostic,
which is consistent with Icepack applying growth as a split
temperature-solve-plus-`thickness_changes` repartitioning step rather than as a
single continuous moving-grid diffusion solve.

The metric-aware thermodynamic solve must use the conservative
`MutableVerticalDiscretization` balance

```math
\mathcal{J}^{n+1} \Delta r_k E_k^{n+1}
= \mathcal{J}^n \Delta r_k E_k^n
+ \delta z_{g,k+1/2} E^{up}_{k+1/2}
- \delta z_{g,k-1/2} E^{up}_{k-1/2}
+ \Delta t \left[(F^E + I)_{k+1/2}^{n+1}
- (F^E + I)_{k-1/2}^{n+1}\right],
```

with the analogous update for prognostic bulk salinity. Here
``\delta z_g = z_g^{n+1} - z_g^n`` comes from the moving Jacobian metric and
``E^{up}`` is the upwind layer value swept by the moving face. Strict BL99
agreement also requires applying CICE's split thickness-change remap after that
metric-aware thermodynamic solve.

## CICE Smoke Gate

The current fixed-thickness local CICE smoke cases are:

```text
/private/tmp/cice-validation/cases_fixed/bl99_cold_conductive_relaxation_MU71_fixed_thickness_balbot
/private/tmp/cice-validation/cases_fixed/bl99_cold_conductive_relaxation_bubbly_fixed_thickness_balbot
```

They were generated with the validation templates in
`validation/cice_bitz_lipscomb/cice_templates`, with
`CICE_BALANCE_COLD_BOTTOM_FLUX=1`, using:

```bash
RUN_CICE_CASES_EXECUTE=1 \
RUN_CICE_CASES_RUN=1 \
RUN_CICE_CASES_BUILD=1 \
CICE_ALLOW_TEMPLATE_OVERWRITE=1 \
CASE_ROOT=/private/tmp/cice-validation/cases_fixed \
CICE_CASE_SUFFIX=_fixed_thickness_balbot \
CICE_INTERNAL_SNOW_THICKNESS=0.00 \
CICE_BALANCE_COLD_BOTTOM_FLUX=1 \
CICE_CASE_FILTER=cold_conductive_relaxation \
validation/cice_bitz_lipscomb/run_cice_cases.sh
```

`CICE_INTERNAL_SNOW_THICKNESS=0.00` applies a scoped local source patch to
CICE's internal initial-condition path so this case starts without snow. The
same runner can restore the default Barrow-like internal snow thickness with
`CICE_INTERNAL_SNOW_THICKNESS=0.20`.

The surface-ablation case additionally uses `calc_Tsfc = .false.` in
`set_nml.bl99_case_surface_ablation` and
`CICE_ABLATION_SURFACE_FLUX_W_M2=75.0` to patch CICE's prescribed `fsurfn_d`
values in `cicecore/cicedyn/general/ice_flux.F90`. The latest forced ablation
case is:

```text
/private/tmp/cice-validation/cases_nosnow_suite/bl99_surface_ablation_MU71_nosnow_ablate_forced
```

with run log:

```text
/private/tmp/cice-validation/cases_nosnow_suite/bl99_surface_ablation_MU71_nosnow_ablate_forced/logs/cice.runlog.260522-110205
```

The controlled surface-warming cases also use `calc_Tsfc = .false.` and patch
prescribed `fsurfn_d` to `8.0 W m^-2`, which warms the interior while limiting
top ablation below the transition threshold. The latest warming cases are:

```text
/private/tmp/cice-validation/cases_nosnow_suite/bl99_surface_warming_MU71_nosnow_warming_forced8
/private/tmp/cice-validation/cases_nosnow_suite/bl99_surface_warming_bubbly_nosnow_warming_forced8
```

The local executable was linked with Homebrew `gfortran`,
`nf-config`, and `nc-config`:

```text
/private/tmp/cice-validation/runs/bl99_cold_conductive_relaxation_MU71_fixed_thickness_balbot/cice
/private/tmp/cice-validation/runs/bl99_cold_conductive_relaxation_bubbly_fixed_thickness_balbot/cice
```

The successful fixed-thickness run logs are:

```text
/private/tmp/cice-validation/cases_fixed/bl99_cold_conductive_relaxation_MU71_fixed_thickness_balbot/logs/cice.runlog.260522-121821
/private/tmp/cice-validation/cases_fixed/bl99_cold_conductive_relaxation_bubbly_fixed_thickness_balbot/logs/cice.runlog.260522-121909
```

Those runs wrote instantaneous hourly, daily, and initial-condition CICE history
files under:

```text
/private/tmp/cice-validation/runs/bl99_cold_conductive_relaxation_MU71_fixed_thickness_balbot/history
/private/tmp/cice-validation/runs/bl99_cold_conductive_relaxation_bubbly_fixed_thickness_balbot/history
```

The parser reads the active central T cell from `iceh_grid.nc`, then extracts
the initial state and 120 hourly records from the CICE history files. The
current fixed-cold per-case smoke summaries are:

```text
validation/cice_bitz_lipscomb/results/cold_conductive_relaxation_MU71_fixed_thickness_cice_smoke_summary.csv
validation/cice_bitz_lipscomb/results/cold_conductive_relaxation_bubbly_fixed_thickness_cice_smoke_summary.csv
validation/cice_bitz_lipscomb/results/cice_case_comparison_summary.csv
```

For the fixed MU71 smoke run it found active cell `j=3, i=3`, initial
`hi = 1.0 m`, `hs = 0 m`, `Tsfc = -20.15 C`, and final day-5
`hi = 1.000003 m`, `hs = 0 m`, `Tsfc = -22.22582 C`. The bubbly fixed run
ends with `hi = 0.9999988 m`, `hs = 0 m`, and `Tsfc = -21.61519 C`.
The summaries also include hourly surface and ocean flux diagnostics such as
`fsurf_ai`, `fcondtop_ai`, `fbot`, `siflcondbot`, `fhocn`, `fswabs`,
`fswthru`, `fsens`, `flat`, `flwdn`, and `flwup`; the first hourly MU71 record
has `fsurf_ai = -48.25298 W m^-2`, `fcondtop_ai = -48.25298 W m^-2`, and
`fhocn = -34.80659 W m^-2`. The first hourly bubbly record has
`fhocn = -37.66438 W m^-2`. CICE's hourly `siflcondbot` output is sparse for
this stream, so the replay reconstructs the BL99 temperature-solve bottom
conductive flux from the parsed state and source formula.

The analyzer also maps each CICE ice-temperature and salinity profile into a
7-layer ClimaSeaIce column state. Its current metrics are:

```text
validation/cice_bitz_lipscomb/results/climaseaice_state_replay_summary.csv
```

The current state-replay gate passes with maximum salinity-profile error
`4.80033825134e-07 ppt`, maximum column-ingest temperature error
`7.1054273576e-15 C`, and fixed-relation temperature round-trip error
`3.5527136788e-15 C`. Reconstructed shifted internal energy from CICE
`Tinz/Sinz` agrees with the ClimaSeaIce relation to `7.45058059692e-09 J m^-3`
in the cold MU71 history.

A deliberately limited fixed-grid forced replay is recorded in:

```text
validation/cice_bitz_lipscomb/results/fixed_grid_forced_replay_summary.csv
```

It uses hourly CICE top fluxes de-area-weighted by `aice`, plus the parsed
CICE bottom ocean temperature as a one-sided conductive boundary. The
fixed-thickness cold run balances the computed bottom conductive flux before `thickness_changes`;
the history has no
reported basal growth and `delta_hi = 3.0e-6 m` over five days. The
reset-each-interval diagnostic starts each interval from the parsed CICE state
and meets the 1% temperature target: maximum reset error is
`0.00608785103382 C`, with maximum relative temperature error
`0.00103279653852`, per-level shifted internal-energy error
`0.00871573484266`, and column-energy error `1.27790080481e-05`.
The free-running fixed-grid replay now also meets the 1% temperature target for
fixed cold with maximum relative temperature error `0.00799365465354`; its
per-level shifted internal-energy diagnostic remains above 1% at
`0.0534201130377`.

The source-level Icepack matrix replay for the same fixed-cold intervals passes
with maximum relative temperature error `7.96658313105e-07`, shifted
internal-energy error `1.1539840022e-06`, and column-energy error
`5.4741495334e-07`.

The analyzer also reports a CICE-thickness-forced replay using the general
ClimaSeaIce fixed-grid step. This diagnostic still uses the CICE heat-flux
history and does not prognose thickness itself, but it conservatively remaps
layer-integrated energy onto each hourly CICE thickness and recomputes the
prescribed salinity profile from normalized depth. For positive basal-growth
intervals it now uses the BL99 basal-congelation fill enthalpy from
`thickness_changes`: CICE computes `qbot` from bottom salinity and `Tbot`, grows
the bottom layer, then repartitions with `adjust_enthalpy`. The validation
replay now uses the core `column_energy_thickness_remap!` helper for this
energy repartitioning and salinity reset. This is distinct
from the `add_new_ice` open-water/frazil distribution path. The source-level
Icepack matrix version of this diagnostic passes the growth cases, so the
remaining general-step growth drift is in ClimaSeaIce's fixed-grid temperature
solve. Ablation still has relative-only and per-level internal-energy drift in
the sequential source-level matrix diagnostic, narrowing that remaining
mismatch to the full BL99 ablation iteration/history-forcing path.

For ablation, conservative remapping brings layer-temperature and column-energy
reset diagnostics below 1% in the ClimaSeaIce fixed-grid replay, but per-level
shifted internal energy remains above the target: `0.0195520787001` for `MU71`
and `0.0196562450921` for `bubbly`. The source-level Icepack matrix plus remap
passes the same ablation internal-energy diagnostic at `0.002322920769` and
`0.0023133988512`. Sequentially applying the same source-level matrix while
following CICE thickness still exceeds the relative-only diagnostic at
`0.0137447491872` and `0.0134220179813`, but passes the official temperature
gate with acceptance ratios `0.877625726166` and `0.866657092284`. The
source-level prognostic-thickness replay now passes ablation after
de-area-weighting CICE `fsurf_ai` and `fcondtop_ai` by `aice`; the general
ClimaSeaIce prognostic ablation path still fails the temperature gate at
`3.57820733413` and `3.49785711177`. The source-level matrix diagnostic uses
Icepack's `nitermax = 100` iteration cap, so this residual is not a
reduced-iteration convergence artifact.
For the MU71 ablation history, the maximum sequential source-level error is
localized to day `5.0`, layer `1` from the top: `-1.29458669567 C` predicted
versus `-1.27703418115 C` in CICE. The final bottom-layer error is
`0.000556099736938 C`, pointing to accumulated top-ablation repartitioning
drift rather than a bulk column energy error. For controlled warming, the
fixed-grid remapped internal-energy
diagnostics pass:
`0.00089974495413` and `0.000915045669336`.

Important limitations of this smoke gate:

- It is a CICE setup/build/run/parse/state-replay gate plus one reset
  fixed-grid forced replay for selected local CICE histories. The strict
  source-level BL99 melting/freezing comparison now passes, while the general
  ClimaSeaIce semi-implicit moving-thickness path remains diagnostic.
- `MutableVerticalDiscretization` metric integration is now covered by the
  column-energy proof gates and reported in the CICE-thickness replay. The
  moving-grid equations are necessary, but CICE-equivalent growth and ablation
  also need CICE-compatible top-ablation and basal-growth enthalpy
  repartitioning.
- The ablation transition is now present and conservative thickness remapping
  brings reset temperature and column energy below 1% for both conductivity
  variants, but per-level shifted internal energy remains above 1%.
- The CICE-thickness-forced replay is diagnostic only: it follows CICE `hi`
  rather than predicting melt/growth from ClimaSeaIce surface and basal
  residuals.
- The `basal_growth` setup intentionally leaves the winter bottom conductive
  flux unbalanced to retain measurable basal congelation; the fixed-cold setup
  uses the validation hook instead.
- It uses idealized `atm_data_type = 'calm'` and `ocn_data_type = 'calm'`.
- This no-snow smoke case uses the runner's scoped `hsno_init = 0.00` source
  patch for CICE's `ice_ic = 'internal'` path; a restart-based initial
  condition would be cleaner for long-term validation.

## Reproduction Commands

Focused local proof gates:

```bash
TEST_GROUP=column_energy julia --project=. test/runtests.jl
```

CICE case setup:

```bash
validation/cice_bitz_lipscomb/run_cice_cases.sh
```

The CICE script defaults to dry-run mode. Set `RUN_CICE_CASES_EXECUTE=1` after
checking machine/compiler settings. On this local macOS/Homebrew environment,
the runner defaults to the validation-provided `conda/homebrew` machine files.

CICE smoke history parser:

```bash
julia --startup-file=no --project=. validation/cice_bitz_lipscomb/analyze_cice_comparison.jl
```

Set `CICE_HISTORY_DIR` to parse a different CICE history directory.
