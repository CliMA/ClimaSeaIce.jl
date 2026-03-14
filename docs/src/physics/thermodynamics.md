# Sea Ice Thermodynamics

Sea ice thermodynamics governs the freezing and melting of ice through heat exchange
with the atmosphere above and the ocean below. ClimaSeaIce implements a **slab sea ice model**,
which represents the ice as a single, vertically-integrated layer with uniform properties.

## The slab model

In the slab approximation, the ice is characterized by its thickness ``h``, concentration ``\aleph``,
and a single temperature at the top surface ``T_u``. The bottom of the ice sits at the
ice-ocean interface, where temperature equals the melting point.

The evolution of ice thickness follows the Stefan condition at both interfaces:
```math
\frac{\partial h}{\partial t} = w_u + w_b
```
where ``w_u`` and ``w_b`` are the interface velocities at the top (atmosphere-ice) and bottom
(ice-ocean) interfaces, respectively.

## Phase transitions and the liquidus

The melting temperature of sea ice depends on its salinity through the **liquidus** relationship.
ClimaSeaIce uses a linear liquidus:
```math
T_m(S) = T_0 - m S
```
where ``T_m`` is the melting temperature, ``T_0`` is the freshwater melting temperature (0°C by default),
``S`` is the salinity, and ``m`` is the liquidus slope (0.054 °C/psu by default).

The [`PhaseTransitions`](@ref) struct encapsulates the thermodynamic properties governing
ice-water transitions:

```@example thermodynamics
using ClimaSeaIce.SeaIceThermodynamics: PhaseTransitions

phase_transitions = PhaseTransitions()
```

### Temperature-dependent latent heat

The latent heat of fusion ``\mathscr{L}`` represents the energy required to transform ice into
liquid water (or released during freezing). It varies with temperature:
```math
\rho_i \mathscr{L}(T) = \rho_i \mathscr{L}_0 + (\rho_\ell c_\ell - \rho_i c_i)(T - T_0)
```
where ``\rho_i`` and ``\rho_\ell`` are the ice and liquid densities, ``c_i`` and ``c_\ell``
are the respective heat capacities, ``\mathscr{L}_0`` is the reference latent heat at
temperature ``T_0``, and ``T`` is the current temperature.

```@example thermodynamics
using ClimaSeaIce.SeaIceThermodynamics: latent_heat

# Latent heat at different temperatures
L_at_0C = latent_heat(phase_transitions, 0.0)
L_at_minus10C = latent_heat(phase_transitions, -10.0)
println("Latent heat at 0°C: ", L_at_0C, " J/m³")
println("Latent heat at -10°C: ", L_at_minus10C, " J/m³")
```

## Heat fluxes and the Stefan condition

The slab model computes thickness changes from the balance of heat fluxes at each interface.

### Internal conductive flux

Within the ice slab, heat is conducted from the warmer bottom to the colder top surface.
The conductive flux is:
```math
Q_i = -k \frac{T_u - T_b}{h}
```
where ``k`` is the thermal conductivity (default 2 W/(m·K) for freshwater ice),
``T_u`` is the top surface temperature, ``T_b`` is the bottom temperature (at the melting point),
and ``h`` is the ice thickness.

The [`ConductiveFlux`](@ref) struct implements this:

```@example thermodynamics
using ClimaSeaIce.SeaIceThermodynamics: ConductiveFlux

# Thermal conductivity of 2 W/(m·K)
conductive_flux = ConductiveFlux(Float64; conductivity = 2.0)
```

### Interface velocities

At each interface, the difference between incoming and outgoing heat fluxes drives
melting or freezing:
```math
w_u = \frac{Q_x - Q_i}{\mathscr{L}(T_u)}, \quad w_b = \frac{Q_i - Q_b}{\mathscr{L}(T_b)}
```
where:
- ``Q_x`` is the external (atmospheric) heat flux into the top surface
- ``Q_i`` is the internal conductive flux
- ``Q_b`` is the oceanic heat flux at the bottom

Negative velocities indicate melting; positive velocities indicate freezing.

## Top surface boundary conditions

The top surface temperature ``T_u`` can be determined in two ways.

### Melting-constrained flux balance

The default approach, [`MeltingConstrainedFluxBalance`](@ref), solves for the temperature
that equilibrates the external and internal heat fluxes:
```math
\sum_n Q_{x,n}(T_u) - Q_i(T_u) = 0
```
subject to the constraint that the surface cannot exceed the melting temperature:
```math
T_u \leq T_m(S)
```

When the surface would otherwise be warmer than the melting point, it is clamped to ``T_m(S)``
and the excess heat drives surface melting.

```@example thermodynamics
using ClimaSeaIce.SeaIceThermodynamics: MeltingConstrainedFluxBalance

top_bc = MeltingConstrainedFluxBalance()
```

### Prescribed temperature

Alternatively, the surface temperature can be prescribed directly using [`PrescribedTemperature`](@ref):

```@example thermodynamics
using ClimaSeaIce.SeaIceThermodynamics: PrescribedTemperature

# Fix the surface at -5°C
top_bc = PrescribedTemperature(-5.0)
```

## Bottom boundary condition

The default bottom boundary condition, `IceWaterThermalEquilibrium`, assumes the ice-ocean
interface is at the salinity-dependent melting temperature:
```math
T_b = T_m(S)
```

This is appropriate when the ocean mixed layer is well-mixed and maintains thermal
equilibrium with the ice bottom.

## External heat fluxes

External heat fluxes drive the thermodynamic evolution of sea ice.

### Radiative emission

The [`RadiativeEmission`](@ref) flux represents longwave radiation emitted by the ice surface:
```math
Q_r = \varepsilon \sigma (T + T_{\text{ref}})^4
```
where ``\varepsilon`` is the emissivity, ``\sigma`` is the Stefan-Boltzmann constant,
and ``T_{\text{ref}}`` converts from Celsius to Kelvin.

### Custom flux functions

The [`FluxFunction`](@ref) wrapper allows arbitrary user-defined heat fluxes that can
depend on position, time, surface temperature, and model fields:

```@example thermodynamics
using ClimaSeaIce.SeaIceThermodynamics: FluxFunction

# A simple sensible heat flux depending on surface temperature
sensible_heat(Tu, clock, fields, parameters) = parameters.coefficient * (parameters.T_air - Tu)
flux = FluxFunction(sensible_heat; parameters = (coefficient = 15.0, T_air = -10.0))
```

## Putting it together: SlabThermodynamics

The [`SlabThermodynamics`](@ref) struct (also aliased as `SlabSeaIceThermodynamics`)
combines all thermodynamic components:

```@example thermodynamics
using Oceananigans
using ClimaSeaIce.SeaIceThermodynamics: SlabThermodynamics, ConductiveFlux, MeltingConstrainedFluxBalance

grid = RectilinearGrid(size=(10, 10), x=(0, 1e5), y=(0, 1e5), topology=(Periodic, Periodic, Flat))

thermodynamics = SlabThermodynamics(grid;
    top_heat_boundary_condition = MeltingConstrainedFluxBalance(),
    internal_heat_flux = ConductiveFlux(Float64; conductivity = 2.0))

thermodynamics
```

When coupled to a [`SeaIceModel`](@ref), these thermodynamics automatically compute thickness
tendencies based on the configured heat fluxes and boundary conditions.

## Snow thermodynamics

ClimaSeaIce supports a **layered snow model** where a snow layer sits on top of the ice slab.
Snow insulates the ice by adding thermal resistance in series, reducing the conductive flux
and slowing ice growth.

### Snow layer physics

When snow is present, the combined conductive flux through the snow+ice slab uses
resistors in series:
```math
F_c = \frac{T_b - T_u}{h_s / k_s + h_i / k_i}
```
where ``h_s, k_s`` are the snow thickness and conductivity, and ``h_i, k_i`` are the ice
thickness and conductivity. Snow typically has ``k_s \approx 0.31`` W/(m·K) compared to
``k_i \approx 2`` W/(m·K) for ice, making it an effective insulating layer.

The snow surface temperature ``T_u`` is determined by balancing the atmospheric heat
fluxes against the combined conductive flux, using the same
[`MeltingConstrainedFluxBalance`](@ref) approach as bare ice. The melting point for snow
is 0°C (rather than the salinity-dependent liquidus for ice).

### Interface temperature

The snow-ice interface temperature ``T_{si}`` is computed analytically from the surface
temperature using the thermal resistance ratio:
```math
T_{si} = T_b + (T_u - T_b) \frac{R_i}{R_s + R_i}
```
where ``R_i = h_i / k_i`` and ``R_s = h_s / k_s``. This temperature is passed to the ice
thermodynamics as a prescribed boundary condition, so the ice layer is completely
unaware of the snow layer above it.

### Snow processes

The snow model includes three processes that modify the snow thickness each time step:

1. **Surface melt**: When the surface temperature reaches the melting point, any excess
   conductive flux first melts snow before reaching the ice. The snow absorbs energy
   ``\rho_s L_s \Delta h_s`` before excess energy melts ice from the top.

2. **Snowfall accumulation**: Snow precipitates at a prescribed rate (kg/m²/s) and
   accumulates only where ice is present (``\aleph > 0``).

3. **Snow-ice formation (flooding)**: When the snow load pushes the ice surface below
   the waterline (negative freeboard), seawater floods the snow layer and converts it
   to ice. The freeboard criterion is:
   ```math
   h_f = h_i (1 - \rho_i / \rho_w) - h_s \rho_s / \rho_w < 0
   ```
   The flooding step is **energy-conserving**: it uses the liquid density ``\rho_w``
   (consistent with the `latent_heat` definition) to compute the mass transfer between
   snow and ice.

### Snow volume conservation

Snow thickness ``h_s`` represents the thickness over the ice-covered fraction ``\aleph``.
The total snow energy per unit area is ``\rho_s L_s h_s \aleph``, consistent with how ice
energy is ``\mathscr{L} h \aleph``. When ice concentration changes (e.g., during
consolidation), the snow thickness is adjusted to conserve the snow volume
``h_s \times \aleph``.

### Setting up a snow-covered model

Use [`SlabSnowThermodynamics`](@ref) for convenient construction with snow defaults:

```@example thermodynamics
using ClimaSeaIce

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

snow_thermodynamics = SlabSnowThermodynamics(grid)

model = SeaIceModel(grid;
    snow_thermodynamics,
    snow_precipitation = 3e-6) # kg/m²/s

set!(model, h=1, ℵ=1, hs=0.1) # 10 cm snow on 1 m ice
```

The [`SeaIceModel`](@ref) constructor automatically wires the layered coupling:
the snow's internal heat flux is replaced with the combined snow+ice conductive flux,
and the ice's top boundary condition becomes a prescribed interface temperature.
