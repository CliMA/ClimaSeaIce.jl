# Layered snow + ice thermodynamics — implementation notes

This page documents the snow-surface energy balance and the closed-form implicit solve for end-of-step ice concentration used by
`_layered_thermodynamic_time_step!`. It assumes familiarity with the slab mass balance described in the Thermodynamics chapter.

## Column geometry

Each cell contains an ice-covered fraction $\aleph$; snow sits on the ice only. The prognostic variables are

- $h_i$ — ice thickness on the ice-covered fraction (m, per-ice)
- $\aleph$ — ice concentration (unitless, per-cell area fraction)
- $h_s$ — snow thickness on the ice-covered fraction (m, per-ice)

with the ice-volume invariant $V = h_i \, \aleph$ tracked per unit *cell* area. Snow volume per unit cell area is $h_s \, \aleph$ and is conserved
across changes in $\aleph$ by the kernel's $h_s \leftarrow h_s \, \aleph^n / \aleph^{n+1}$ rescale.

## Per-ice vs. per-cell fluxes

Two conventions coexist and must be treated consistently:

- $Q_{ui}$ (top external heat flux) and $Q_{bi}$ (bottom external heat flux) are delivered by the coupler **per unit cell area** (radiation and
  turbulent fluxes $\times \aleph$ on the ice path, interface heat $\times \aleph$ on the ocean path).
- $Q_{is}$ (column conductive flux $(T_b - T_u)/R$, with $R = h_s/k_s + h_i/k_i$) is intrinsically **per unit ice area**: the thermal resistance only
  applies where ice exists.
- $Q_{ii}$ (ice-only internal conductive flux) is similarly per-ice.

## Snow-surface energy balance

The snow surface sits on the ice-covered fraction only, so the surface energy balance is a **per-ice** balance:

```math
\delta Q = \frac{Q_{ui}}{\aleph^n} - Q_{is}, \qquad
Q_s = \min\!\left(\max(0,\, -\delta Q),\; \frac{\rho_s \, \mathscr{L} \, h_s^n}{\Delta t}\right).
```

$Q_s$ (positive when melting) is the per-ice latent power absorbed by snow melt, giving a thickness rate $G_s^- = Q_s / (\rho_s \, \mathscr{L})$.

The atmospheric flux $Q_{ui}$ passed by the coupler is per-cell, so we divide by $\aleph^n$ before comparing with $Q_{is}$.
When $\aleph^n = 1$ this reduces to the original formulation; when $\aleph^n < 1$ it avoids the spurious factor of $\aleph^n$ that
would otherwise suppress the snow melt rate.

## Implicit concentration update

The slab mass balance and the `ProportionalEvolution` concentration rule are both linear in the end-of-step concentration $\aleph^{n+1}$:

```math
\partial_t V = \frac{(Q_{ui} + Q_s \, \aleph^{n+1}) - Q_{bi}}{\rho_i \, \mathscr{L}} \;\equiv\; \alpha + \beta \, \aleph^{n+1},
```

with

```math
\alpha = \frac{Q_{ui} - Q_{bi}}{\rho_i \, \mathscr{L}}, \qquad
\beta  = \frac{Q_s}{\rho_i \, \mathscr{L}},
```

and

```math
\aleph^{n+1} = \aleph^n + \Delta t \, C \, \partial_t V,
```

with

```math
C = \begin{cases}
\dfrac{\aleph^n}{2 \, h_i^n} & \text{if } \partial_t V < 0 \quad \text{(melt)} \\[6pt]
\dfrac{1 - \aleph^n}{h^c}    & \text{if } \partial_t V \ge 0 \quad \text{(freeze)}.
\end{cases}
```

Defining $K = \Delta t \, C$, the fixed point $\aleph^{n+1} \, (1 - K \beta) = \aleph^n + K \alpha$ has the **closed-form solution**

```math
\boxed{\,\aleph^{n+1} = \frac{\aleph^n + K \alpha}{1 - K \beta}\,}.
```

## Branch selection

The correct branch (melt vs. freeze) depends on the sign of $\partial_t V$ at the solution, which isn't known up front. The kernel computes *both*
branches, $\aleph^m$ (melt) and $\aleph^f$ (freeze), and picks the one where $\operatorname{sign}(\partial_t V(\aleph^{n+1}))$ is consistent with the branch assumption:

- if $\alpha + \beta \, \aleph^m < 0$, use the melt-branch solution $\aleph^m$
- otherwise, use the freeze-branch solution $\aleph^f$

This is branchless (via `ifelse`) and correct at the $\alpha \approx 0$ boundary where the first-guess sign could be misleading. On ocean scales
$|K \beta| \lesssim 10^{-6}$, so the sign rarely flips between guess and solution; the test is defensive for the atypical case where atmospheric heating
almost exactly cancels the column conductive flux.

## NaN guards

All divisions are protected:

- $Q_{ui} / \aleph^n$ falls back to $0$ if $\aleph^n \le 0$
- $\aleph^n / (2 \, h_i^n)$ falls back to $0$ if $h_i^n \le 0$
- $(1 - \aleph^n) / h^c$ falls back to $0$ if $h^c \le 0$
- $1 / (1 - K \beta)$ falls back to the numerator $\aleph^n + K \alpha$ if the denominator is within $\varepsilon$ of zero

The $Q_s$ cap $\rho_s \, \mathscr{L} \, h_s^n / \Delta t$ keeps $K \beta$ comfortably below 1 in practice; the guard is a defensive fallback rather than a routine path.

## Edge cases deferred to `ice_volume_update`

The closed-form solves the implicit $\aleph^{n+1}$ update but doesn't apply volume clipping, ridging, or pathological-case handling. After the
solve, the kernel calls `ice_melt_freeze_tendency` once at the self-consistent $Q_{ui}^{\mathrm{eff}} = Q_{ui} + Q_s \, \aleph^{n+1}$ and passes the resulting
$\partial_t V$ through `ice_volume_update`, which handles:

- $V^{n+1} = \max(0,\, V^n + \Delta t \, \partial_t V)$ when the ice melts completely
- ridging cap $\aleph^+ \to 1$ with compensating thickness adjustment
- degenerate states ($\aleph \le 0$, $\partial_t V = 0$, $h^+ = 0$)

## Energy conservation

The implicit solve eliminates the per-step $Q_s \, (1 - \aleph^{n+1}) \, \Delta t \, A$ leak that an explicit or single-pass formulation would incur.
