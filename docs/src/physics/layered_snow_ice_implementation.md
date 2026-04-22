# Layered snow + ice thermodynamics — implementation notes

This page documents the snow-surface energy balance and the closed-form implicit solve for end-of-step ice concentration used by
`_layered_thermodynamic_time_step!`. It assumes familiarity with the slab mass balance described in the Thermodynamics chapter.

## Column geometry

Each cell contains an ice-covered fraction `ℵ`; snow sits on the ice only. The prognostic variables are

- `h`  — ice thickness on the ice-covered fraction (m, per-ice)
- `ℵ`  — ice concentration (unitless, per-cell area fraction)
- `hs` — snow thickness on the ice-covered fraction (m, per-ice)

with the ice-volume invariant `V = h · ℵ` tracked per unit *cell* area. Snow volume per unit cell area is `hs · ℵ` and is conserved 
across changes in `ℵ` by the kernel's `hs ← hs · ℵⁿ/ℵⁿ⁺¹` rescale.

## Per-ice vs. per-cell fluxes

Two conventions coexist and must be treated consistently:

- `Qui` (top external heat flux), `Qbi` (bottom external heat flux) are delivered by the coupler **per unit cell area** (radiation and
  turbulent fluxes × ℵ on the ice path, interface heat × ℵ on the ocean path).
- `Qis` (column conductive flux `(Tb − Tu)/R`, `R = hs/ks + hi/ki`) is intrinsically **per unit ice area**: the thermal resistance only
  applies where ice exists.
- `Qii` (ice-only internal conductive flux) is similarly per-ice.

## Snow-surface energy balance

The snow surface sits on the ice-covered fraction only, so the surface energy balance is a **per-ice** balance:

```math
\delta Q = \frac{Qui}{\aleph^n} - Qis, \qquad
Q_s = \min\left(\max(0, -\delta Q), \frac{\rho_s \mathcal{L}_s h_s^n}{\Delta t}\right).
```

`Qs` (positive when melting) is the per-ice latent power absorbed by snow melt, giving a thickness rate `Gs⁻ = Qs/(ρₛ ℒₛ)`.

The atmospheric flux `Qui` passed by the coupler is per-cell, so we divide by `ℵⁿ` before comparing with `Qis`. 
When `ℵⁿ = 1` this reduces to the original formulation; when `ℵⁿ < 1` it avoids the spurious factor of `ℵⁿ` that 
would otherwise suppress the snow melt rate.

## Implicit concentration update

The slab mass balance and the `ProportionalEvolution` concentration rule are both linear in the end-of-step concentration `ℵⁿ⁺¹`:

```math
\partial_t V = \frac{(Qui + Q_s \aleph^{n+1}) - Qbi}{\rho_i \mathcal{L}}  \equiv \alpha + \beta\, \aleph^{n+1},
```

with

```math
\alpha = \frac{Qui - Qbi}{\rho_i \mathcal{L}}, \quad
\beta  = \frac{Q_s}{\rho_i \mathcal{L}},
```

and

```math
\aleph^{n+1} = \aleph^n + \Delta t\, C\, \partial_t V,
```

with

```math
C = \begin{cases}
\frac{\aleph^n}{2 h^n} & \text{if } \partial_t V < 0 \ \text{(melt)} \\
\frac{1 - \aleph^n}{h^c} & \text{if } \partial_t V \ge 0 \ \text{(freeze)}.
\end{cases}
```

Defining `K = Δt · C`, the fixed point `ℵⁿ⁺¹ · (1 − K β) = ℵⁿ + K α` has the **closed-form solution**

```math
\boxed{\aleph^{n+1} = \frac{\aleph^n + K \alpha}{1 - K \beta}}.
```

## Branch selection

The correct branch (melt vs. freeze) depends on the sign of `∂t_V` at the solution, which isn't known up front. The kernel computes *both*
branches and picks the one where `sign(∂t_V(ℵⁿ⁺¹))` is consistent with the branch assumption:

- if `α + β · ℵᵐ < 0`, use the melt-branch solution `ℵᵐ`
- otherwise, use the freeze-branch solution `ℵᶠ`

This is branchless (via `ifelse`) and correct at the `α ≈ 0` boundary where the first-guess sign could be misleading. On ocean scales
`|K β| ≲ 10⁻⁶`, so the sign rarely flips between guess and solution; the test is defensive for the atypical case where atmospheric heating
almost exactly cancels the column conductive flux.

## NaN guards

All divisions are protected:

- `Qui / ℵⁿ`      → falls back to `0` if `ℵⁿ ≤ 0`
- `ℵⁿ / (2 hⁿ)`   → falls back to `0` if `hⁿ ≤ 0`
- `(1 − ℵⁿ) / hᶜ` → falls back to `0` if `hᶜ ≤ 0`
- `1 / (1 − K β)` → falls back to the numerator `ℵⁿ + K α` if the denominator is within `eps` of zero

The `Qs` cap `ρₛ ℒₛ hₛⁿ / Δt` keeps `K β` comfortably below 1 in practice; the guard is a defensive fallback rather than a routine path.

## Edge cases deferred to `ice_volume_update`

The closed-form solves the implicit ℵⁿ⁺¹ update but doesn't apply volume clipping, ridging, or pathological-case handling. After the
solve, the kernel calls `ice_melt_freeze_tendency` once at the self-consistent `Quie = Qui + Qs · ℵⁿ⁺¹` and passes the resulting
`∂t_V` through `ice_volume_update`, which handles:

- `V^{n+1} = max(0, V^n + Δt ∂t_V)` when the ice melts completely
- ridging cap `ℵ⁺ → 1` with compensating thickness adjustment
- degenerate states (`ℵ ≤ 0`, `∂t_V = 0`, `h⁺ = 0`)

## Energy conservation

The implicit solve eliminates the per-step `Qs · (1 − ℵⁿ⁺¹) · Δt · A`. leak that an explicit or single-pass formulation would incur. 
