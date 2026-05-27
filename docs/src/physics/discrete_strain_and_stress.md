# Metric-consistent strain rates and stress divergence

The EVP rheology closes a loop between velocity and internal force through two operators:

1. the **strain rate** ``\dot{\boldsymbol\varepsilon}``, which feeds the stress ``\boldsymbol\sigma``;
2. the **stress divergence** ``\nabla\cdot\boldsymbol\sigma``, which feeds the velocity tendency.

This page discretizes the strain rate, then *derives* the discrete stress divergence by
demanding that it be the discrete adjoint of the strain operator. The pay-off is a single
algebraic identity — discrete summation by parts — from which the metric factors of
``\nabla\cdot\boldsymbol\sigma`` fall out automatically, and which guarantees that the internal
stress can only remove kinetic energy on *any* grid.

## Variable placement

Velocities sit on an Arakawa C-grid; stresses are split between cell centers and corners.

| Quantity | Symbol | Location |
|---|---|---|
| zonal velocity | ``u`` | `(Face, Center)` |
| meridional velocity | ``v`` | `(Center, Face)` |
| normal stresses | ``\sigma_{11}, \sigma_{22}`` | `(Center, Center)` |
| shear stress | ``\sigma_{12}`` | `(Face, Face)` |

The deformation invariants follow the same placement: divergence and tension at centers, shear
at corners. Grid spacings are the curvilinear scale factors; we write ``\Delta x_\ell``,
``\Delta y_\ell`` for the spacing at location ``\ell\in\{ccc, fcc, cfc, ffc\}`` and
``A_\ell = \Delta x_\ell\,\Delta y_\ell`` for the cell area.

## The one tool: discrete summation by parts

The difference operators shift a field between staggered locations,
```math
(\delta_i^{f\to c} p)_i = p_{i+1} - p_i ,\qquad (\delta_i^{c\to f} a)_i = a_i - a_{i-1},
```
(in code, `δxᶜᵃᵃ` and `δxᶠᵃᵃ`; analogously ``\delta_j`` with `δyᵃᶜᵃ`, `δyᵃᶠᵃ`). Reindexing a
periodic sum, or one in which the fields vanish at the boundary, gives

```math
\boxed{\;\sum_i (\delta_i^{f\to c} p)_i\, a_i = -\sum_i p_i\, (\delta_i^{c\to f} a)_i\;}
\tag{SBP}
```

i.e. the two staggered differences are negative adjoints. This is the only identity used below.

## Step 1 — strain rates, from the continuum

Work in orthogonal curvilinear coordinates ``(q_1, q_2)`` — the grid directions ``i, j`` — with
scale factors ``\Delta x, \Delta y`` (the physical length per unit coordinate, so that
``\mathrm ds^2 = \Delta x^2\,\mathrm dq_1^2 + \Delta y^2\,\mathrm dq_2^2``; on the grid these are
the cell spacings, which may vary in space). For physical velocity components ``(u, v)`` the
physical components of the rate-of-strain tensor
``\dot\varepsilon_{kl} = \tfrac12(\partial_k u_l + \partial_l u_k)`` are

```math
\dot\varepsilon_{11} = \frac{1}{\Delta x}\,\partial_i u + \frac{v}{\Delta x\,\Delta y}\,\partial_j \Delta x,
\qquad
\dot\varepsilon_{22} = \frac{1}{\Delta y}\,\partial_j v + \frac{u}{\Delta x\,\Delta y}\,\partial_i \Delta y,
```
```math
\dot\varepsilon_{12} = \frac12\left[\frac{\Delta y}{\Delta x}\,\partial_i\!\Big(\frac{v}{\Delta y}\Big)
                                 + \frac{\Delta x}{\Delta y}\,\partial_j\!\Big(\frac{u}{\Delta x}\Big)\right].
```

The terms in ``\partial \Delta x, \partial \Delta y`` are the metric (curvature) terms; they
vanish only when the scale factors are constant. Form the three invariants
``\dot\varepsilon_D = \dot\varepsilon_{11} + \dot\varepsilon_{22}`` (divergence),
``\dot\varepsilon_T = \dot\varepsilon_{11} - \dot\varepsilon_{22}`` (tension), and
``\dot\varepsilon_S = 2\dot\varepsilon_{12}`` (shear). The product-rule identities

```math
\frac{\partial_i(\Delta y\, u)}{\Delta x\,\Delta y} = \frac{1}{\Delta x}\partial_i u + \frac{u}{\Delta x\,\Delta y}\partial_i\Delta y,
\qquad
\frac{\Delta y^2\,\partial_i(u/\Delta y)}{\Delta x\,\Delta y} = \frac{1}{\Delta x}\partial_i u - \frac{u}{\Delta x\,\Delta y}\partial_i\Delta y,
```

(and their ``\partial_j`` counterparts) fold each curvature term into one weighted derivative.
Adding and subtracting then collapses the invariants to

```math
\dot\varepsilon_D = \frac{1}{\Delta x\,\Delta y}\Big[\partial_i(\Delta y\,u) + \partial_j(\Delta x\,v)\Big]
\quad(\text{divergence}),
```
```math
\dot\varepsilon_T = \frac{1}{\Delta x\,\Delta y}\Big[\Delta y^2\,\partial_i(u/\Delta y) - \Delta x^2\,\partial_j(v/\Delta x)\Big]
\quad(\text{tension}),
```
```math
\dot\varepsilon_S = \frac{1}{\Delta x\,\Delta y}\Big[\Delta x^2\,\partial_j(u/\Delta x) + \Delta y^2\,\partial_i(v/\Delta y)\Big]
\quad(\text{shear}).
```

The ``\Delta x^2, \Delta y^2`` weights are therefore not a modeling choice: they *are* the
curvature terms of the strain tensor, repackaged by the identities above. For constant scale
factors they cancel and the invariants reduce to ``\partial_x u + \partial_y v``,
``\partial_x u - \partial_y v`` and ``\partial_y u + \partial_x v``.

Discretize by replacing ``\partial`` with the staggered difference that lands on the target
location, and evaluating each scale factor where its field lives. Divergence and tension land
at centers, shear at corners:

```math
\dot\varepsilon_D = \frac{1}{A_{ccc}}\Big[\delta_i(\Delta y_{fcc}\,u) + \delta_j(\Delta x_{cfc}\,v)\Big],
```
```math
\dot\varepsilon_T = \frac{1}{A_{ccc}}\Big[\Delta y_{ccc}^2\,\delta_i(u/\Delta y_{fcc}) - \Delta x_{ccc}^2\,\delta_j(v/\Delta x_{cfc})\Big],
```
```math
\dot\varepsilon_S = \frac{1}{A_{ffc}}\Big[\Delta x_{ffc}^2\,\delta_j(u/\Delta x_{fcc}) + \Delta y_{ffc}^2\,\delta_i(v/\Delta y_{cfc})\Big].
```

The tensor components are recovered pointwise,
``\dot\varepsilon_{11}=\tfrac12(\dot\varepsilon_D+\dot\varepsilon_T)``,
``\dot\varepsilon_{22}=\tfrac12(\dot\varepsilon_D-\dot\varepsilon_T)``,
``\dot\varepsilon_{12}=\tfrac12\dot\varepsilon_S``, and the EVP scalar is
``\Delta=\sqrt{\dot\varepsilon_D^2 + e^{-2}(\dot\varepsilon_T^2+\dot\varepsilon_S^2)}``.

## Step 2 — the discrete internal-stress power

Define the area-weighted inner product at a location ``\ell``,
``\langle a,b\rangle_\ell = \sum_\ell A_\ell\, a\, b``. The continuous internal-stress power is
``W=\int \boldsymbol\sigma:\dot{\boldsymbol\varepsilon}\,\mathrm dA``. Writing
``\sigma_I=\sigma_{11}+\sigma_{22}``, ``\sigma_{II}=\sigma_{11}-\sigma_{22}`` and using
``\boldsymbol\sigma:\dot{\boldsymbol\varepsilon}= \tfrac12\sigma_I\dot\varepsilon_D + \tfrac12\sigma_{II}\dot\varepsilon_T + \sigma_{12}\dot\varepsilon_S``,
the discrete power placed at the natural locations of each invariant is

```math
W = \tfrac12\langle \sigma_I,\dot\varepsilon_D\rangle_{ccc}
  + \tfrac12\langle \sigma_{II},\dot\varepsilon_T\rangle_{ccc}
  + \langle \sigma_{12},\dot\varepsilon_S\rangle_{ffc}.
\tag{W}
```

## Step 3 — derive the divergence by duality

We *define* the discrete ``\nabla\cdot\boldsymbol\sigma`` to be the operator for which the power
delivered to the velocity field equals ``-W``:

```math
\langle (\nabla\cdot\boldsymbol\sigma)_x,\, u\rangle_{fcc}
+ \langle (\nabla\cdot\boldsymbol\sigma)_y,\, v\rangle_{cfc}
\;=\; -\,W .
\tag{D}
```

The right-hand side is linear in ``(u,v)``; collecting the coefficient of each ``u`` and ``v``
gives the operator. Take the three ``u``-dependent pieces of ``-W`` and move every difference off
``u`` with (SBP). The area weights cancel against the ``1/A_{ccc}``, ``1/A_{ffc}`` in the strain
rates:

**Divergence piece.**
```math
-\tfrac12\langle\sigma_I,\dot\varepsilon_D\rangle_{ccc}\big|_u
= -\tfrac12\sum_{ccc}\sigma_I\,\delta_i(\Delta y_{fcc}u)
\;\overset{\text{(SBP)}}{=}\; +\tfrac12\sum_{fcc} (\Delta y_{fcc}u)\,\delta_i\sigma_I
= \big\langle \tfrac{1}{A_{fcc}}\cdot\tfrac12\,\Delta y_{fcc}\,\delta_i\sigma_I,\; u\big\rangle_{fcc}.
```

**Tension piece.**
```math
-\tfrac12\sum_{ccc}\sigma_{II}\,\Delta y_{ccc}^2\,\delta_i(u/\Delta y_{fcc})
\;\overset{\text{(SBP)}}{=}\;
+\tfrac12\sum_{fcc}\frac{u}{\Delta y_{fcc}}\,\delta_i(\Delta y_{ccc}^2\sigma_{II})
= \Big\langle \tfrac{1}{A_{fcc}}\cdot\tfrac12\,\frac{\delta_i(\Delta y_{ccc}^2\sigma_{II})}{\Delta y_{fcc}},\;u\Big\rangle_{fcc}.
```

**Shear piece.**
```math
-\sum_{ffc}\sigma_{12}\,\Delta x_{ffc}^2\,\delta_j(u/\Delta x_{fcc})
\;\overset{\text{(SBP)}}{=}\;
+\sum_{fcc}\frac{u}{\Delta x_{fcc}}\,\delta_j(\Delta x_{ffc}^2\sigma_{12})
= \Big\langle \tfrac{1}{A_{fcc}}\cdot\frac{\delta_j(\Delta x_{ffc}^2\sigma_{12})}{\Delta x_{fcc}},\;u\Big\rangle_{fcc}.
```

Reading off the coefficient of ``u`` (and, identically, of ``v``):

```math
\boxed{\;(\nabla\cdot\boldsymbol\sigma)_x = \frac{1}{A_{fcc}}\!\left[
\tfrac12\,\Delta y_{fcc}\,\delta_i\sigma_I
+ \tfrac12\,\frac{\delta_i(\Delta y_{ccc}^2\,\sigma_{II})}{\Delta y_{fcc}}
+ \frac{\delta_j(\Delta x_{ffc}^2\,\sigma_{12})}{\Delta x_{fcc}}\right]\;}
```
```math
\boxed{\;(\nabla\cdot\boldsymbol\sigma)_y = \frac{1}{A_{cfc}}\!\left[
\tfrac12\,\Delta x_{cfc}\,\delta_j\sigma_I
- \tfrac12\,\frac{\delta_j(\Delta x_{ccc}^2\,\sigma_{II})}{\Delta x_{cfc}}
+ \frac{\delta_i(\Delta y_{ffc}^2\,\sigma_{12})}{\Delta y_{cfc}}\right]\;}
```

Three things were *derived*, not chosen:

- the squared scale factors ``\Delta y_{ccc}^2``, ``\Delta x_{ffc}^2``, … in the divergence are
  forced by the same factors in the strain rate, through (SBP);
- the ``\tfrac12`` on the ``\sigma_I,\sigma_{II}`` terms and the coefficient ``1`` on the
  ``\sigma_{12}`` term come straight from the contraction weights in (W);
- the **minus** sign on the ``\sigma_{II}`` term of ``(\nabla\cdot\boldsymbol\sigma)_y`` is the
  minus sign of ``\dot\varepsilon_T``'s ``v``-part (``\sigma_{22}=(\sigma_I-\sigma_{II})/2``).

## The property they satisfy

By construction (D), the discrete strain operator ``\mathsf E:(u,v)\mapsto(\dot\varepsilon_D,\dot\varepsilon_T,\dot\varepsilon_S)``
and the discrete divergence are a negative-adjoint pair,
```math
\langle \nabla\cdot\boldsymbol\sigma,\, \boldsymbol u\rangle
= -\,\big[\tfrac12\langle\sigma_I,\dot\varepsilon_D\rangle + \tfrac12\langle\sigma_{II},\dot\varepsilon_T\rangle + \langle\sigma_{12},\dot\varepsilon_S\rangle\big]
= -\langle \boldsymbol\sigma, \dot{\boldsymbol\varepsilon}\rangle ,
```
for arbitrary ``\boldsymbol\sigma`` and ``\boldsymbol u``. This is a discrete integration by parts
with **no truncation error and no dependence on grid uniformity** — the scale factors and areas
are arbitrary positive numbers throughout the derivation.

The consequence for the momentum equation ``m\,\partial_t\boldsymbol u = \nabla\cdot\boldsymbol\sigma + \dots``
is exact discrete control of the internal-stress energy budget:
```math
\frac{\mathrm d}{\mathrm d t}\sum \tfrac12 m\,|\boldsymbol u|^2 A \;\Big|_{\text{internal}}
= \langle \nabla\cdot\boldsymbol\sigma,\,\boldsymbol u\rangle
= -\langle \boldsymbol\sigma, \dot{\boldsymbol\varepsilon}\rangle \;\le\; 0,
```
the last inequality because the visco-plastic stress satisfies
``\boldsymbol\sigma:\dot{\boldsymbol\varepsilon}\ge 0``. The internal stress is therefore
*strictly dissipative at the discrete level*: it cannot inject kinetic energy into grid-scale
modes. A pair of operators that does not satisfy (D) loses this guarantee on a non-uniform grid,
where the missing ``\Delta x^2/\Delta y^2`` weights leave a residual of indefinite sign.

### Numerical verification

Identity (D) is *exact*, so it is easy to check directly. Build the strain and divergence
operators on a curvilinear grid (where ``\Delta x`` genuinely varies in space — a
`LatitudeLongitudeGrid` will do), place arbitrary fields on the velocity and stress points, and
compare the two scalars. Taper the fields to zero away from the latitude boundaries so the
boundary terms in (SBP) vanish; then the only thing being tested is the operator pair
(`divergence`/`tension`/`shear` and `∂ⱼ_σ₁ⱼ`/`∂ⱼ_σ₂ⱼ` are the functions from the table below,
written here in standalone form that takes the stress fields directly):

```julia
# u, v on (Face,Center)/(Center,Face); σ₁₁, σ₂₂ on centers; σ₁₂ on corners — all windowed.
LHS = 0.0; RHS = 0.0
for j in 5:Ny-4, i in 1:Nx                       # interior; periodic in i
    LHS += Azᶠᶜᶜ(i,j,1,g)*u[i,j,1]*∂ⱼ_σ₁ⱼ(i,j,1,g, s11,s22,s12)
    LHS += Azᶜᶠᶜ(i,j,1,g)*v[i,j,1]*∂ⱼ_σ₂ⱼ(i,j,1,g, s11,s22,s12)
    RHS -= Azᶜᶜᶜ(i,j,1,g)*(0.5*σI(i,j,1,g,s11,s22) *divergence(i,j,1,g,u,v)
                         + 0.5*σII(i,j,1,g,s11,s22)*tension(i,j,1,g,u,v))
    RHS -= Azᶠᶠᶜ(i,j,1,g)* s12[i,j,1]            *shear(i,j,1,g,u,v)
end
@show LHS, RHS, abs(LHS - RHS) / max(abs(LHS), abs(RHS))
```

On a latitude–longitude grid this returns a relative mismatch of order ``10^{-15}`` — the two
sides agree to round-off, confirming the operators are an adjoint pair *independently of the
metric*. The same check on a uniform Cartesian grid would pass for the naive operators too; it is
precisely the curvilinear case that distinguishes them.

## Mapping to the code

| Object | Function | File |
|---|---|---|
| ``\dot\varepsilon_D,\dot\varepsilon_T,\dot\varepsilon_S`` | `ice_divergence`, `ice_tension`, `ice_shear` | `src/Rheologies/elasto_visco_plastic_rheology.jl` |
| tensor components | `strain_rate_xx`, `_yy`, `_xy` | same |
| ``\zeta`` at center and corner | `_compute_evp_viscosities!` | same |
| ``\sigma_I,\sigma_{II}`` | `σI`, `σII` | `src/Rheologies/ice_stress_divergence.jl` |
| ``(\nabla\cdot\sigma)_x,(\nabla\cdot\sigma)_y`` | `∂ⱼ_σ₁ⱼ`, `∂ⱼ_σ₂ⱼ` | same |

The ``q/\Delta`` and ``\Delta^2 q`` building blocks appear in the code as the helpers
`q_over_Δyᶠᶜᶜ`, `Δy²_qᶜᶜᶜ`, etc., which are exactly the per-location scale-factor weights that
make the strain and divergence operators an adjoint pair.
