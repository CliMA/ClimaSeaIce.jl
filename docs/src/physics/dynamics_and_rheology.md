# Sea Ice Dynamics and Rheology

Sea ice dynamics describes the horizontal motion of ice in response to wind stress,
ocean currents, Coriolis forces, and internal stresses arising from ice-ice interactions.
The **rheology** determines how sea ice deforms and transmits stress—a crucial ingredient
for capturing the characteristic fractures, leads, and ridges observed in real sea ice.

## The momentum equation

The sea ice momentum equation describes the balance of forces acting on an ice element:
```math
\frac{\partial \boldsymbol{u}}{\partial t} + \boldsymbol{f} \times \boldsymbol{u} =
    \frac{1}{m_i} \boldsymbol{\nabla} \cdot \boldsymbol{\sigma} +
    \frac{\boldsymbol{\tau}_a}{m_i} + \frac{\boldsymbol{\tau}_o}{m_i}
```
where:
- ``\boldsymbol{u} = (u, v)`` is the ice velocity
- ``\boldsymbol{f}`` is the Coriolis parameter
- ``m_i = \rho_i h \aleph`` is the ice mass per unit area (density × thickness × concentration)
- ``\boldsymbol{\sigma}`` is the internal stress tensor
- ``\boldsymbol{\tau}_a`` is the atmospheric stress (wind drag)
- ``\boldsymbol{\tau}_o`` is the oceanic stress (water drag)

The divergence of the stress tensor, ``\boldsymbol{\nabla} \cdot \boldsymbol{\sigma}``, represents
internal forces arising from ice deformation.

## Setting up ice dynamics

The [`SeaIceMomentumEquation`](@ref) struct configures the dynamical components:

```@example dynamics
using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics: SeaIceMomentumEquation, SemiImplicitStress

L = 512kilometers
grid = RectilinearGrid(size = (64, 64),
                       x = (0, L),
                       y = (0, L),
                       topology = (Bounded, Bounded, Flat))

# Ocean velocity field for ice-ocean drag
Uₒ = XFaceField(grid)
Vₒ = YFaceField(grid)

# Cyclonic pattern; ∂u/∂y + ∂v/∂x = 0
𝓋ₒ = 0.01 # m/s
set!(Uₒ, (x, y) -> 𝓋ₒ * (2y - L) / L)
set!(Vₒ, (x, y) -> 𝓋ₒ * (L - 2x) / L)

dynamics = SeaIceMomentumEquation(grid;
    coriolis = FPlane(f = 1e-4),
    bottom_momentum_stress = SemiImplicitStress(uₑ = Uₒ, vₑ = Vₒ))

dynamics
```

## External stresses

### Semi-implicit stress formulation

The ice experiences drag from both the atmosphere above and the ocean below.
The [`SemiImplicitStress`](@ref) implements a quadratic drag law:
```math
\tau_u = \rho_e C_D \sqrt{(u_e - u_i^n)^2 + (v_e - v_i^n)^2} (u_e - u_i^{n+1})
```
```math
\tau_v = \rho_e C_D \sqrt{(u_e - u_i^n)^2 + (v_e - v_i^n)^2} (v_e - v_i^{n+1})
```
where ``\rho_e`` is the external fluid density, ``C_D`` is the drag coefficient,
``(u_e, v_e)`` are the external (ocean or atmosphere) velocities, and
``(u_i, v_i)`` are the ice velocities.

The "semi-implicit" treatment uses velocities from the current time step (``n``) to compute
the magnitude, while the velocity difference uses the new time step (``n+1``), improving
numerical stability.

```@example dynamics
using ClimaSeaIce.SeaIceDynamics: SemiImplicitStress

# Ocean drag (default parameters for seawater)
ocean_stress = SemiImplicitStress(
    uₑ = Uₒ,
    vₑ = Vₒ,
    ρₑ = 1026.0,  # seawater density [kg/m³]
    Cᴰ = 5.5e-3   # drag coefficient
)
```

For atmospheric stress, use lower density and appropriate drag coefficient:
```@example dynamics
Uₐ = XFaceField(grid)
Vₐ = YFaceField(grid)

atm_stress = SemiImplicitStress(
    uₑ = Uₐ,
    vₑ = Vₐ,
    ρₑ = 1.3,     # air density [kg/m³]
    Cᴰ = 1.2e-3   # air-ice drag coefficient
)
```

## Rheology: how ice deforms

The rheology determines the relationship between ice deformation (strain rates) and
internal stresses. Sea ice exhibits complex mechanical behavior: it can flow viscously
under slow deformation, fracture plastically under large stresses, and respond
elastically on short timescales.

### The stress tensor

The internal stress tensor ``\boldsymbol{\sigma}`` is a function of the strain rate tensor
``\dot{\boldsymbol{\varepsilon}}``:
```math
\sigma_{ij} = \sigma_{ij}(\dot{\varepsilon}_{kl})
```

The strain rates are computed from velocity gradients:
```math
\dot{\varepsilon}_{11} = \frac{\partial u}{\partial x}, \quad
\dot{\varepsilon}_{22} = \frac{\partial v}{\partial y}, \quad
\dot{\varepsilon}_{12} = \frac{1}{2}\left(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}\right)
```

The specific form of ``\sigma_{ij}(\dot{\varepsilon}_{kl})`` depends on the chosen rheology.

### Viscous rheology

The simplest rheology is purely viscous, where stress is linearly proportional to strain rate:
```math
\sigma_{ij} = 2\nu \dot{\varepsilon}_{ij}
```
with constant viscosity ``\nu``.

```@example dynamics
using ClimaSeaIce.Rheologies: ViscousRheology

# Viscous rheology with ν = 1000 m²/s
rheology = ViscousRheology(ν = 1000.0)
```

This is computationally simple but doesn't capture the plastic behavior that creates
realistic ice features like leads and ridges.

### Elasto-Visco-Plastic (EVP) rheology

The [`ElastoViscoPlasticRheology`](@ref) is the workhorse for realistic sea ice simulations.

Physically, sea ice behaves as a **visco-plastic** material: it creeps viscously under
low stress and yields plastically when stress exceeds a threshold. However, solving the
visco-plastic equations directly requires implicit nonlinear solvers, which are computationally
expensive.

The EVP method introduces a small amount of **artificial elasticity** as a numerical device.
This elasticity allows the stress equations to be subcycled explicitly within each time step,
iterating toward the visco-plastic solution. When sufficiently converged, the elastic terms
vanish and the solution satisfies the true visco-plastic rheology.

#### The constitutive relation

The EVP rheology uses the following stress-strain relationship:
```math
\sigma_{ij} = 2\eta \dot{\varepsilon}_{ij} + \left[(\zeta - \eta)(\dot{\varepsilon}_{11} + \dot{\varepsilon}_{22}) - \frac{P}{2}\right] \delta_{ij}
```
where:
- ``\dot{\varepsilon}_{ij}`` is the strain rate tensor
- ``\eta`` is the shear viscosity
- ``\zeta`` is the bulk viscosity
- ``P`` is the ice strength (pressure)
- ``\delta_{ij}`` is the Kronecker delta

The shear stress (``\sigma_{12}``) depends only on the shear viscosity ``\eta``, while the
normal stresses (``\sigma_{11}``, ``\sigma_{22}``) depend on both viscosities and the ice strength.

#### Ice strength

The ice strength ``P`` determines when plastic yielding occurs:
```math
P = P_\star h \exp\left[-C(1 - \aleph)\right]
```
where:
- ``P_\star`` is the compressive strength parameter (default: 27,500 N/m²)
- ``h`` is the ice thickness
- ``C`` is the compaction hardening exponent (default: 20)
- ``\aleph`` is the ice concentration

Thick, compact ice is strong; thin, sparse ice is weak.

#### The yield curve

Stresses cannot exceed the **yield curve**, an ellipse in stress space with eccentricity ``e``:
```math
\frac{\sigma_I^2}{P^2} + \frac{\sigma_{II}^2}{(P/e)^2} = 1
```
where ``\sigma_I`` and ``\sigma_{II}`` are principal stresses. The default eccentricity
``e = 2`` means ice yields more easily under shear than compression.

#### Viscosities

The bulk and shear viscosities are computed from the ice strength and a
visco-plastic parameter ``\Delta``:
```math
\zeta = \frac{P}{2\Delta}, \quad \eta = \frac{\zeta}{e^2}
```
where:
```math
\Delta = \sqrt{\delta^2 + e^{-2} s^2}
```
with ``\delta = \dot{\varepsilon}_{11} + \dot{\varepsilon}_{22}`` (divergence) and
``s = \sqrt{(\dot{\varepsilon}_{11} - \dot{\varepsilon}_{22})^2 + 4\dot{\varepsilon}_{12}^2}`` (shear).

A minimum ``\Delta_{\min}`` prevents infinite viscosities in quiescent regions.

```@example dynamics
using ClimaSeaIce.Rheologies: ElastoViscoPlasticRheology

rheology = ElastoViscoPlasticRheology(
    ice_compressive_strength = 27500,  # P* [N/m²]
    ice_compaction_hardening = 20,     # C
    yield_curve_eccentricity = 2,      # e
    minimum_plastic_stress = 2e-9      # Δ_min
)
```

#### Dynamic substepping

The EVP formulation uses dynamic substepping following Kimmritz et al. (2016).
Rather than fixed substeps, a spatially-varying relaxation parameter ``\alpha`` accelerates
convergence where the ice is weak and allows slower relaxation where it is strong:

```math
\sigma_{ij}^{p+1} = \sigma_{ij}^{p} + \frac{\sigma_{ij}^{p+1} - \sigma_{ij}^{p}}{\alpha}
```

The relaxation parameter is computed from a stability criterion:
```math
\alpha = \text{clamp}\left(\sqrt{\frac{\zeta c_\alpha \Delta t}{m_i A_z}}, \alpha^-, \alpha^+\right)
```
where ``c_\alpha`` controls the relaxation strength and ``\alpha^-``, ``\alpha^+`` are
the minimum and maximum bounds.

```@example dynamics
using ClimaSeaIce.Rheologies: ElastoViscoPlasticRheology

# Full EVP configuration
rheology = ElastoViscoPlasticRheology(
    ice_compressive_strength = 27500,
    ice_compaction_hardening = 20,
    yield_curve_eccentricity = 2,
    minimum_plastic_stress = 2e-9,
    min_relaxation_parameter = 50,   # α⁻
    max_relaxation_parameter = 300,  # α⁺
    relaxation_strength = π^2        # c_α
)
```

## Momentum solvers

### Split-explicit solver

The default [`SplitExplicitSolver`](@ref) subcycles the momentum equation within each
thermodynamic time step. This is necessary because the EVP rheology requires many
iterations to converge:

```@example dynamics
using ClimaSeaIce.SeaIceDynamics: SplitExplicitSolver

# 150 momentum substeps per time step
solver = SplitExplicitSolver(grid; substeps = 150)
```

The solver iterates:
1. Compute strain rates from current velocities
2. Update viscosities and ice strength
3. Compute stresses using the EVP relaxation
4. Update velocities from the momentum equation

### Explicit solver

For simpler rheologies (like viscous), a single explicit update per time step suffices:

```@example dynamics
using ClimaSeaIce.SeaIceDynamics: ExplicitSolver

solver = ExplicitSolver()
```

## Complete example

Here's a complete setup for ice dynamics with EVP rheology:

```@example dynamics
using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics: SeaIceMomentumEquation, SemiImplicitStress

Lx = Ly = 512kilometers
grid = RectilinearGrid(size = (64, 64),
                       x = (0, Lx),
                       y = (0, Ly),
                       topology = (Bounded, Bounded, Flat))

# Ocean drag
Uₒ = XFaceField(grid)
Vₒ = YFaceField(grid)
set!(Uₒ, (x, y) -> 0.01 * (2y - 512e3) / 512e3)
set!(Vₒ, (x, y) -> 0.01 * (512e3 - 2x) / 512e3)

# Atmospheric stress fields (to be updated during simulation)
τₐu = XFaceField(grid)
τₐv = YFaceField(grid)

dynamics = SeaIceMomentumEquation(grid;
    coriolis = FPlane(f = 1e-4),
    rheology = ElastoViscoPlasticRheology(),
    solver = SplitExplicitSolver(grid; substeps = 150),
    top_momentum_stress = (u = τₐu, v = τₐv),
    bottom_momentum_stress = SemiImplicitStress(uₑ = Uₒ, vₑ = Vₒ)
)

model = SeaIceModel(grid;
    dynamics = dynamics,
    ice_thermodynamics = nothing,  # Pure dynamics
    advection = WENO(order = 7))

set!(model, h = 0.3, ℵ = 1.0)

model
```

This configuration produces realistic ice deformation patterns including
linear kinematic features (LKFs) that form when ice fractures under stress.
