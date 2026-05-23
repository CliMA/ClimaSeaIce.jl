# Bitz-Lipscomb CICE Validation

This page tracks validation of the fixed-salinity column-energy
thermodynamics against CICE/Icepack BL99 column physics.

The current status is source matching, pre-CICE proof gates, and local no-snow
CICE column runs with parsed instantaneous history output, diagnostic
ClimaSeaIce state replay, and reset-each-interval fixed-grid forced replay. The
column energy solve now has a `MutableVerticalDiscretization` path that evolves
energy and prognostic salinity with the moving vertical Jacobian and explicit
moving-face enthalpy/salinity flux. The CICE comparison scripts now report both
the older fixed-grid/remap diagnostics and a moving-metric replay. The strict
BL99 validation gate now passes for all eight no-snow histories by using the
source-level Icepack BL99 temperature matrix, CICE-style thickness-change
remapping, and the source-level prognostic-thickness replay. The current
`cold_conductive_relaxation` CICE runs use a validation-only
balanced-bottom-flux hook and pass the fixed-thickness transition gate for both
`MU71` and `bubbly`; the fixed-thickness replay now holds the predicted
thickness fixed while still reporting CICE's observed `hi` drift separately.
Matching the Icepack BL99 `cp_ice = 2106` and `cp_ocn = 4218` constants brings
fixed-cold, basal-growth, and controlled surface-warming reset temperature,
per-level shifted internal energy, and column-energy diagnostics below 1%.
A ClimaSeaIce source-level Icepack BL99 temperature-matrix replay followed by
CICE-style conservative thickness remapping now passes all thickness-changing
reset diagnostics, including forced surface ablation. The forced
surface-ablation cases produce the intended top melt for both conductivity
variants; their raw fixed-grid reset replays miss the 1% target, and
ClimaSeaIce fixed-grid remapping still leaves per-level shifted internal energy
above 1% for ablation. The general ClimaSeaIce CICE-thickness-forced replay
passes fixed-cold, controlled warming, and basal growth by the official
temperature gate, but forced ablation still fails in the semi-implicit
split-solve path. The moving-metric replay confirms that the conservative
`MutableVerticalDiscretization` Jacobian equations are necessary, but CICE
agreement for growth and ablation still comes from Icepack's split
`thickness_changes` enthalpy redistribution rather than from a single
continuous moving-grid diffusion solve.

## Reference

The inspected local checkout is `/private/tmp/cice-validation/CICE`.

| component | commit |
| --- | --- |
| CICE | `3c2064659ccf693e644e0cfd21140565ddb8e355` |
| Icepack | `68097c59aceb2040f8c0cef189bccfe7590a8a32` |

## Matched Physics

The CICE/Icepack salinity profile is implemented as

```math
S(s) = \frac{S_{\max}}{2}
\left[1 - \cos\left(\pi s^{a / (s + b)}\right)\right],
```

where `s` is normalized downward depth from the top ice surface,
`S_max = 3.2`, `a = 0.407`, and `b = 0.573`.

The fixed-salinity brine-pocket relation uses the shifted internal energy

```math
E(T, S) =
\rho_i c_i T
- \rho_i (c_w - c_i) \mu S
- \rho_i L_0 \frac{\mu S}{T}.
```

The CICE-compatible conductivity closures are:

```math
K_i^{MU71} =
\max\left(K_0 + \beta \frac{S}{\min(-\epsilon, T)}, K_{min}\right),
```

and

```math
K_i^{bubbly} =
\max\left[
\frac{\rho_i}{917}
\left(2.11 - 0.011 T + 0.09 \frac{S}{\min(-\epsilon, T)}\right),
K_{min}
\right].
```

CICE stores interface `kh` in conductance units. For equal layer thickness
`h`, CICE's interior ice value is

```math
k_h =
\frac{2 K_{k-1} K_k}{(K_{k-1} + K_k) h}.
```

ClimaSeaIce stores face conductivity for the grid-operator form, so the
corresponding equal-spacing face conductivity is the harmonic mean

```math
K_f =
\frac{2 K_{k-1} K_k}{K_{k-1} + K_k}.
```

For positive basal growth, the diagnostic remap follows Icepack's BL99
`thickness_changes` convention: basal congelation uses bottom salinity and
`Tbot` to define the bottom growth enthalpy, grows the bottom layer, then
repartitions energy onto equal layers with the `adjust_enthalpy` overlap remap.
That path is separate from CICE's open-water/frazil `add_new_ice` distribution.
For positive top ablation, Icepack `thickness_changes` removes ice using
`qm = zqin` for `ktherm != 2`, then repartitions the remaining layer enthalpy
onto equal layers with the same `adjust_enthalpy` overlap remap.

The public ClimaSeaIce names used for this BL99 validation are
`FixedDrainedIceSalinityProfile`, `FixedSalinityBrinePocketEnergyRelation`,
`MaykutUntersteinerConductivity`, `BubblyBrineConductivity`,
`MeltingConstrainedSurfaceFluxBalance`, and
`OceanFreezingTemperatureBoundary`. The last two are descriptive aliases for
the existing `MeltingConstrainedFluxBalance` and `IceWaterThermalEquilibrium`
boundary conditions.

## Current Dashboard

```@raw html
<table>
<tr><th>Gate</th><th>Status</th></tr>
<tr><td>Fixed salinity profile</td><td>pass in <code>TEST_GROUP=column_energy</code></td></tr>
<tr><td>Fixed salinity relation</td><td>pass in <code>TEST_GROUP=column_energy</code></td></tr>
<tr><td>Conductivity closures</td><td>pass in <code>TEST_GROUP=column_energy</code></td></tr>
<tr><td>Moving-grid column metric</td><td>pass in <code>TEST_GROUP=column_energy</code>; no-flux mutable-grid steps conserve layer-integrated energy and salinity</td></tr>
<tr><td>CICE cold conductive relaxation</td><td>fixed-thickness transition and reset replay pass for <code>conduct='MU71'</code> and <code>conduct='bubbly'</code></td></tr>
<tr><td>CICE history parser</td><td>pass for parsed central-column histories</td></tr>
<tr><td>ClimaSeaIce state replay</td><td>pass for profile mapping and relation round trip</td></tr>
<tr><td>Icepack temperature-matrix replay</td><td>source-level replay plus CICE-style conservative remap passes all reset histories</td></tr>
<tr><td>Reset fixed-grid forced replay</td><td>temperature and column energy pass for fixed-cold, basal-growth, and controlled warming histories; raw forced ablation fails</td></tr>
<tr><td>Conservative thickness remap</td><td>ClimaSeaIce fixed-grid remap passes temperature and column energy; per-level internal energy still fails for ablation diagnostics</td></tr>
<tr><td>Icepack matrix with CICE thickness</td><td>sequential source-level replay passes the official temperature gate for all histories; forced ablation still exceeds the relative-only 1% diagnostic and per-level internal-energy target</td></tr>
<tr><td>CICE-thickness-forced replay</td><td>passes fixed-cold, controlled warming, and basal growth by temperature; general ClimaSeaIce path still fails forced ablation</td></tr>
<tr><td>Moving-metric CICE-thickness replay</td><td>passes fixed-cold by temperature and controlled warming by temperature/internal energy; growth and ablation still need CICE-compatible boundary-motion repartitioning</td></tr>
<tr><td>Free-running fixed-grid replay</td><td>diagnostic; fixed-cold temperature passes, moving-thickness cases still require mutable-grid metric coupling</td></tr>
<tr><td>CICE surface warming</td><td>run for <code>MU71</code> and <code>bubbly</code>; controlled warming transition and reset replay pass</td></tr>
<tr><td>CICE surface ablation</td><td>run for <code>MU71</code> and <code>bubbly</code>; top-melt transition and remapped reset replay pass</td></tr>
<tr><td>CICE basal growth</td><td>run for <code>MU71</code> and <code>bubbly</code>; transition and reset replay pass</td></tr>
<tr><td>Strict BL99 validation</td><td>pass for all eight histories through the source-level BL99 matrix prognostic replay and source-level CICE-thickness temperature, enthalpy, and column-energy gates</td></tr>
</table>
```

The source extraction details and scripts live under
`validation/cice_bitz_lipscomb/`.

The current fixed-thickness CICE smoke logs are
`/private/tmp/cice-validation/cases_fixed/bl99_cold_conductive_relaxation_MU71_fixed_thickness_balbot/logs/cice.runlog.260522-121821`
and
`/private/tmp/cice-validation/cases_fixed/bl99_cold_conductive_relaxation_bubbly_fixed_thickness_balbot/logs/cice.runlog.260522-121909`.
The fixed-cold parsed summaries cover the active central T cell `j=3, i=3`,
with the initial state plus 120 hourly records through day 5. The first hourly
MU71 record has `fsurf_ai = -48.25298 W m^-2`,
`fcondtop_ai = -48.25298 W m^-2`, and `fhocn = -34.80659 W m^-2`, where
`fhocn` is the validation hook's balancing bottom heat flux.

The current multi-case CICE aggregate is
`validation/cice_bitz_lipscomb/results/cice_case_comparison_summary.csv`.
The one-figure visual summary is generated by
`validation/cice_bitz_lipscomb/plot_cice_single_column_validation.jl` and
written to
`validation/cice_bitz_lipscomb/results/cice_single_column_validation.png`.

```@raw html
<table>
<tr><th>Case</th><th>Conduct</th><th>Transition</th><th>Raw reset</th><th>Icepack-matrix remapped reset</th><th>Thickness-remapped reset</th><th>CICE-thickness replay</th><th>Note</th></tr>
<tr><td><code>cold_conductive_relaxation</code></td><td><code>MU71</code></td><td>pass</td><td>pass, <code>0.00103279653852</code></td><td>pass, <code>5.14638749576e-06</code></td><td>pass, <code>0.00103304883639</code></td><td>pass, <code>0.00799516820487</code></td><td>Balanced-bottom fixed-thickness run; <code>delta_hi = 3.0e-6 m</code>, <code>delta(hi/aice) = 3.30000099003e-6 m</code>, no basal growth.</td></tr>
<tr><td><code>cold_conductive_relaxation</code></td><td><code>bubbly</code></td><td>pass</td><td>pass, <code>0.00110474853209</code></td><td>pass, <code>5.74253411496e-07</code></td><td>pass, <code>0.00110506393058</code></td><td>pass, <code>0.00778523135678</code></td><td>Balanced-bottom fixed-thickness run; <code>delta_hi = -1.2e-6 m</code>, <code>delta(hi/aice) = -4.00000319978e-7 m</code>, no basal growth.</td></tr>
<tr><td><code>basal_growth</code></td><td><code>MU71</code></td><td>pass</td><td>pass, <code>0.00263944227201</code></td><td>pass, <code>9.26538101245e-06</code></td><td>pass, <code>0.00107403724222</code></td><td>pass, <code>0.00875578804266</code></td><td>Uses the unbalanced winter bottom-flux setup to retain measurable basal congelation.</td></tr>
<tr><td><code>basal_growth</code></td><td><code>bubbly</code></td><td>pass</td><td>pass, <code>0.00283466651198</code></td><td>pass, <code>1.06792979317e-05</code></td><td>pass, <code>0.00113305358948</code></td><td>pass, <code>0.00852913566164</code></td><td>Same basal-growth setup with CICE/Icepack <code>conduct='bubbly'</code>.</td></tr>
<tr><td><code>surface_warming</code></td><td><code>MU71</code></td><td>pass</td><td>pass, <code>0.0012871994515</code></td><td>pass, <code>3.84029189641e-06</code></td><td>pass, <code>0.000357230763411</code></td><td>pass, <code>0.0022686519764</code></td><td>Prescribed <code>8.0 W m^-2</code> surface flux warms the column while limiting ice loss to <code>0.0045599 m</code>.</td></tr>
<tr><td><code>surface_warming</code></td><td><code>bubbly</code></td><td>pass</td><td>pass, <code>0.0012921956637</code></td><td>pass, <code>3.89723767928e-06</code></td><td>pass, <code>0.000363279184033</code></td><td>pass, <code>0.00221754423547</code></td><td>Same controlled warming setup with CICE/Icepack <code>conduct='bubbly'</code>.</td></tr>
<tr><td><code>surface_ablation</code></td><td><code>MU71</code></td><td>pass</td><td>fail, <code>0.0190483763055</code></td><td>pass, <code>0.000919646651979</code></td><td>pass, <code>0.0077718269457</code></td><td>fail, <code>0.0542395682928</code></td><td>Prescribed positive surface flux removes <code>0.10343 m</code> grid-cell ice volume; <code>hi/aice</code> decreases by <code>0.0531144458597 m</code>.</td></tr>
<tr><td><code>surface_ablation</code></td><td><code>bubbly</code></td><td>pass</td><td>fail, <code>0.0191031535505</code></td><td>pass, <code>0.000917650982306</code></td><td>pass, <code>0.007828341593</code></td><td>fail, <code>0.0523392114103</code></td><td>Same forced ablation setup; <code>hi/aice</code> decreases by <code>0.0531160845177 m</code>.</td></tr>
</table>
```

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
all eight histories; the largest
residual is about `3.24e-4 m` in the
forced surface-ablation cases. The integrated top melt is `0.103046700238 m`
for `surface_ablation` with `MU71` and `0.103035348451 m` with `bubbly`; the
integrated basal growth is `0.0518581154167 m` for `basal_growth` with `MU71`
and `0.056307245 m` with `bubbly`.

The aggregate now also reports a first prognostic-thickness forced replay. This
diagnostic applies the CICE-compatible top flux split, routes
the de-area-weighted top-flux split into surface melt, computes basal growth
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

The aggregate also reports `icepack_matrix_prognostic_thickness_forced_*`, a
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

The aggregate also includes shifted internal-energy, CICE negative-enthalpy,
and column-energy error columns. The shifted internal-energy relative floor is
fixed from the cold conductive MU71 case at `2.37722677947e6 J m^-3`. The CICE
negative-enthalpy floor is fixed from the same case at
`2.95867899346e8 J m^-3`, using ``q_i = E - \rho_i L_0``. With the shifted
floor, remapped reset temperature and column energy pass for the fixed-cold,
growth, warming, and ablation reset diagnostics. Per-level shifted internal
energy passes for fixed-cold (`0.00871786496474` for `MU71`,
`0.00929842173781` for `bubbly`), basal growth (`0.00912957581408` for `MU71`,
`0.00963171963647` for `bubbly`), and controlled warming
(`0.00089974495413` and `0.000915045669336`). Remapped surface ablation remains
above 1% in shifted-energy coordinates (`0.0195520787001` and
`0.0196562450921`), but passes as CICE negative enthalpy (`0.000152277126163`
and `0.000153054745318`).

Temperature acceptance ratios use the validation gate directly: relative
temperature error divided by `0.01`, except layers with `|T_CICE| < 2 C` may
instead pass by absolute error divided by `0.02 C`. The aggregate reports
paired `*_temperature_acceptance_status` and
`*_temperature_acceptance_ratio` columns for the reset, remapped reset,
CICE-thickness-forced, moving-metric, source-level matrix, top-ablation reset,
and free-running replay paths.

The aggregate's `icepack_temperature_matrix_reset_*` columns replay the
source-level BL99 temperature matrix with CICE's top conductive flux and bottom
ocean temperature from each interval. For the fixed-cold histories, that replay
matches CICE with maximum relative temperature errors `7.96658313105e-07` for
`MU71` and `1.01110189246e-06` for `bubbly`; the corresponding shifted
internal-energy errors are `1.1539840022e-06` and `2.32760063197e-06`.

The aggregate also reports
`icepack_temperature_matrix_thickness_remapped_reset_*`, which applies the
source-level temperature solve and then the conservative enthalpy remap used by
Icepack `thickness_changes`. This source-level remapped diagnostic passes all
eight histories. For forced ablation, maximum relative temperature errors are
`0.000919646651979` for `MU71` and `0.000917650982306` for `bubbly`; shifted
internal-energy errors are `0.002322920769` and `0.0023133988512`.

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
single continuous moving-grid diffusion solve. The local moving-grid proof gate
therefore includes boundary fill enthalpy support, while the strict CICE replay
uses the source-level split repartitioning sequence.

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

The diagnostic ClimaSeaIce state-replay metrics are in
`validation/cice_bitz_lipscomb/results/climaseaice_state_replay_summary.csv`.
They verify that CICE `Sinz` matches `FixedDrainedIceSalinityProfile` to
`4.80033825134e-07 ppt`, and that the mapped ClimaSeaIce column state round
trips through `FixedSalinityBrinePocketEnergyRelation` to
`3.5527136788e-15 C`. Reconstructed shifted internal energy from CICE
`Tinz/Sinz` agrees with the ClimaSeaIce relation to `7.45058059692e-09 J m^-3`
in the cold MU71 history.

The fixed-grid forced replay metrics are in
`validation/cice_bitz_lipscomb/results/fixed_grid_forced_replay_summary.csv`.
This diagnostic uses hourly CICE top fluxes de-area-weighted by `aice`, plus
the parsed CICE bottom ocean temperature as a one-sided conductive boundary. The
fixed-thickness cold run uses
`CICE_BALANCE_COLD_BOTTOM_FLUX=1`, which preserves the BL99 temperature solve
and feeds the computed bottom conductive flux to `thickness_changes` as the
balancing ocean heat flux. The history has no reported basal growth and
`delta_hi = 3.0e-6 m` over five days. The reset-each-interval replay starts
each interval from the parsed CICE state and has maximum temperature error
`0.00608785103382 C`, maximum relative temperature error `0.00103279653852`,
per-level shifted internal-energy error `0.00871573484266`, and column-energy
error `1.27790080481e-05`. The free-running fixed-grid variant now also meets
the 1% temperature target for fixed cold with maximum relative temperature error
`0.00799365465354`; its per-level shifted internal-energy diagnostic remains
above 1% at `0.0534201130377`.

The source-level Icepack matrix replay for the same fixed-cold intervals passes
with maximum relative temperature error `7.96658313105e-07`, shifted
internal-energy error `1.1539840022e-06`, and column-energy error
`5.4741495334e-07`.

The forced ablation case uses `calc_Tsfc = .false.` and a local CICE source
patch setting prescribed `fsurfn_d` to `75.0 W m^-2`. It passes the transition
gate with about `0.10343 m` ice loss over five days. The raw fixed-grid reset
replay has maximum relative temperature error `0.0190483763055` for `MU71` and
`0.0191031535505` for `bubbly`. Conservatively remapping layer-integrated
energy from the start thickness to each target CICE thickness reduces those
temperature errors to `0.0077718269457` and `0.007828341593`, and column
energy errors to `1.98104938567e-06` and `2.01807627982e-06`. Per-level
shifted internal-energy errors remain above the 1% target:
`0.0195520787001` and `0.0196562450921`. The source-level Icepack matrix plus
conservative remap passes the same ablation internal-energy diagnostic at
`0.002322920769` and `0.0023133988512`. Sequentially applying the same
source-level matrix while following CICE thickness still exceeds the
relative-only diagnostic at `0.0137447491872` and `0.0134220179813`, but passes
the official temperature gate with acceptance ratios `0.877625726166` and
`0.866657092284`. The source-level prognostic-thickness replay now passes
ablation after de-area-weighting CICE `fsurf_ai` and `fcondtop_ai` by `aice`;
the general ClimaSeaIce prognostic ablation path still fails the temperature
gate at `3.57820733413` and `3.49785711177`. This diagnostic uses Icepack's
`nitermax = 100` temperature-matrix iteration cap.
For the MU71 ablation history, the maximum sequential source-level error is
localized to day `5.0`, layer `1` from the top: `-1.29458669567 C` predicted
versus `-1.27703418115 C` in CICE. The final bottom-layer error is
`0.000556099736938 C`, confirming that the drift is concentrated at the
ablating surface.

The controlled warming case uses the same prescribed-flux path with `fsurfn_d`
set to `8.0 W m^-2` for two days. It warms the mean ice temperature by about
`0.011 C`, limits ice loss to `0.0045599 m`, and passes reset replay for both
conductivity variants. The remapped reset maximum relative temperature errors
are `0.000357230763411` for `MU71` and `0.000363279184033` for `bubbly`;
the corresponding internal-energy errors are `0.00089974495413` and
`0.000915045669336`.

This is a setup/build/run/parse/state-replay gate plus selected reset and
thickness-forced diagnostics. The CICE-thickness-forced replays follow the CICE
`hi` sequence instead of predicting melt/growth from ClimaSeaIce surface and
basal residuals. Growth intervals in the remap diagnostic use the BL99
basal-congelation fill enthalpy from Icepack `thickness_changes`; the
validation replay now uses ClimaSeaIce's core
`column_energy_thickness_remap!` helper to repartition energy and reset fixed
salinity from normalized depth. The source-level matrix version passes growth,
so the remaining general-step growth drift is in ClimaSeaIce's temperature
solve and moving-boundary treatment.
Ablation still has relative-only and per-level internal-energy drift even in
the sequential source-level CICE-thickness diagnostic, but it passes the
official temperature gate with the near-melt absolute fallback and has small
CICE-negative-enthalpy and column-energy errors. The strict source-level BL99
validation now passes all eight histories. The general ClimaSeaIce
semi-implicit path remains a diagnostic target for future work: it needs the
prognostic `MutableVerticalDiscretization` Jacobian/face-flux integration plus
CICE-compatible top-ablation and basal-growth enthalpy redistribution to become
CICE-equivalent without invoking the source-level Icepack matrix sequence.
