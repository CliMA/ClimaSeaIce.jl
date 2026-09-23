# ClimaSeaIce.jl

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16143708.svg?style=flat-square)](https://doi.org/10.5281/zenodo.16143708)
[![Aqua](https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![codecov](https://codecov.io/gh/CliMA/ClimaSeaIce.jl/graph/badge.svg?token=3Smw4jVzZG)](https://codecov.io/gh/CliMA/ClimaSeaIce.jl)
[![Documentation](https://img.shields.io/badge/documentation-in%20development-orange)](https://clima.github.io/ClimaSeaIceDocumentation/dev)

ClimaSeaIce.jl is Julia software for simulating the freezing, melting, and horizontal motion of sea ice on CPUs and GPUs. It is designed for climate-scale sea-ice modeling, supports standalone simulations, and can be coupled to ocean models built with [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl).

## What ClimaSeaIce provides

- Sea-ice thermodynamics for freezing, melting, conductive heat transfer, and surface and basal heat-flux parameterizations.
- Sea-ice dynamics with explicit and split-explicit momentum solvers.
- Multiple rheological closures, including viscous and elasto-visco-plastic options.
- Support for standalone sea-ice models and models coupled to Oceananigans-based ocean simulations.

## Installation

Install the latest registered release with

```julia
using Pkg
Pkg.add("ClimaSeaIce")
```

To use the development version from GitHub, use

```julia
using Pkg
Pkg.add(url = "https://github.com/CliMA/ClimaSeaIce.jl")
```

## Documentation

Documentation is available at:

- [Stable documentation](https://clima.github.io/ClimaSeaIceDocumentation/stable/)
- [Development documentation](https://clima.github.io/ClimaSeaIceDocumentation/dev/)

The documentation includes model setup, physics descriptions, examples, and internal implementation notes.

## Quick start

The main entry point is `SeaIceModel`. Thermodynamic parameterizations live under `SeaIceThermodynamics`, and momentum closures and rheologies live under `SeaIceDynamics` and `Rheologies`.

Example scripts are available in [`examples/`](https://github.com/CliMA/ClimaSeaIce.jl/tree/main/examples).

## Citing

If you use ClimaSeaIce.jl in research, teaching, or derived software, please cite the Zenodo record:

> Silvestri, S. et al. (2026). *CliMA/ClimaSeaIce.jl: v0.5.6* (v0.5.6). Zenodo. <https://doi.org/10.5281/zenodo.16143708>
