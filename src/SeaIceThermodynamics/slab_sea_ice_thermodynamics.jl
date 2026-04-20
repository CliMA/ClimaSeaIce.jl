import Oceananigans: fields, prognostic_fields, prognostic_state, restore_prognostic_state!

struct ProportionalEvolution end

struct SlabThermodynamics{ST, HBC, CF, CE}
    top_surface_temperature   :: ST
    heat_boundary_conditions  :: HBC
    internal_heat_flux        :: CF
    concentration_evolution   :: CE
end

Adapt.adapt_structure(to, t::SlabThermodynamics) =
    SlabThermodynamics(Adapt.adapt(to, t.top_surface_temperature),
                       Adapt.adapt(to, t.heat_boundary_conditions),
                       Adapt.adapt(to, t.internal_heat_flux),
                       Adapt.adapt(to, t.concentration_evolution))

const SSIT = SlabThermodynamics

"""
    snow_slab_thermodynamics(grid; kw...)

Construct a `SlabThermodynamics` with default parameters appropriate for snow:
conductivity = 0.31 W/(m K). Bulk density and all phase-transition parameters
live on `SeaIceModel` (as `snow_density` and `phase_transitions` respectively).
"""
function snow_slab_thermodynamics(grid;
                                  conductivity = 0.31,
                                  kw...)

    FT = eltype(grid)
    internal_heat_flux = ConductiveFlux(FT, conductivity = conductivity)
    return SlabThermodynamics(grid; internal_heat_flux, kw...)
end

Base.summary(therm::SSIT) = "SlabThermodynamics"

function Base.show(io::IO, therm::SSIT)
    print(io, "SlabThermodynamics", '\n')
    print(io, "└── top_surface_temperature: ", summary(therm.top_surface_temperature))
end

fields(therm::SSIT) = (; Tu = therm.top_surface_temperature)
prognostic_fields(therm::SSIT) = NamedTuple()

"""
    SlabThermodynamics(grid; kw...)

A minimal slab representation of a single sea-ice or snow layer. Stores
the top surface temperature, top/bottom heat boundary conditions, the raw
internal-flux coefficient (e.g. a `ConductiveFlux`), and the concentration
evolution rule. Phase-transition parameters (densities, latent heats,
liquidus) are stored at the `SeaIceModel` level and threaded into the
tendency kernels via `model.phase_transitions`.
"""
function SlabThermodynamics(grid;
                            top_surface_temperature        = nothing,
                            top_heat_boundary_condition    = MeltingConstrainedFluxBalance(),
                            bottom_heat_boundary_condition = IceWaterThermalEquilibrium(),
                            # Default internal flux: thermal conductivity of 2 kg m s⁻³ K⁻¹, appropriate for freshwater ice
                            internal_heat_flux             = ConductiveFlux(eltype(grid), conductivity=2),
                            concentration_evolution        = ProportionalEvolution())

    if isnothing(top_surface_temperature)
        if top_heat_boundary_condition isa PrescribedTemperature
            top_surface_temperature = top_heat_boundary_condition.temperature
            top_surface_temperature = field((Center, Center, Nothing), top_surface_temperature, grid)
        else
            top_surface_temperature = Field{Center, Center, Nothing}(grid)
        end
    end

    heat_boundary_conditions = (top = top_heat_boundary_condition,
                                bottom = bottom_heat_boundary_condition)

    return SlabThermodynamics(top_surface_temperature,
                              heat_boundary_conditions,
                              internal_heat_flux,
                              concentration_evolution)
end

"""
    sea_ice_slab_thermodynamics(grid; kw...)

Construct a `SlabThermodynamics` with default parameters appropriate for sea ice:
conductivity = 2 W/(m K). Bulk density and all phase-transition parameters live
on `SeaIceModel` (as `sea_ice_density` and `phase_transitions` respectively).
"""
sea_ice_slab_thermodynamics(grid; kw...) = SlabThermodynamics(grid; kw...)

#####
##### Flux-function assembly (used by tendency kernels)
#####

"""
    internal_flux_function(flux, liquidus, bottom_heat_boundary_condition)

Wrap a raw internal-flux coefficient (`ConductiveFlux`, `IceSnowConductiveFlux`,
user `Function`, or user struct) in the `FluxFunction` shape expected by the
surface-temperature solver and by `getflux`. The wrapper is built at the
tendency-kernel level so that `liquidus` and `bottom_heat_boundary_condition`
can be read from `model.phase_transitions` and the slab's heat BCs without
threading those values through `SlabThermodynamics` at construction time.

If `flux` is already a `FluxFunction`, it is returned unchanged — the user
is assumed to have fully assembled the wrapper themselves, and the kernel
does not inject its own parameters.

To plug a custom flux-coefficient struct into the bare-ice slab, extend the
`flux_kernel` dispatch:

```julia
struct MyInternalFlux{T}
    some_parameter :: T
end

# Signature must match the standard FluxFunction kernel:
#   (i, j, grid, Tu, clock, fields, parameters) -> Q
@inline function my_internal_flux(i, j, grid, Tu, clock, fields, parameters)
    flux = parameters.flux         # ::MyInternalFlux
    # ...use `flux.some_parameter`, `fields.h`, etc...
end

ClimaSeaIce.SeaIceThermodynamics.flux_kernel(::MyInternalFlux) = my_internal_flux
```

Once `flux_kernel` dispatches, `MyInternalFlux` can be passed to any
`SlabThermodynamics`/`sea_ice_slab_thermodynamics`/`snow_slab_thermodynamics`
constructor as `internal_heat_flux = MyInternalFlux(...)`.

The layered snow+ice path additionally assumes that each layer's flux
coefficient carries a `.conductivity` field and that the coupling is
resistors-in-series. Custom flux types for a layered column are out of
scope for this refactor.
"""
@inline function internal_flux_function(flux, liquidus, bottom_heat_boundary_condition)
    parameters = (flux = flux,
                  liquidus = liquidus,
                  bottom_heat_boundary_condition = bottom_heat_boundary_condition)

    return FluxFunction(flux_kernel(flux);
                        parameters,
                        top_temperature_dependent = true)
end

# Pass-through when the user has already assembled the `FluxFunction` wrapper.
@inline internal_flux_function(f::FluxFunction, liquidus, bottom_heat_boundary_condition) = f

"""
    flux_kernel(flux)

Return the kernel function that computes the slab's internal heat flux given
a raw flux-coefficient `flux`. This is a public extension point: users who
define a custom flux struct should add a method

```julia
ClimaSeaIce.SeaIceThermodynamics.flux_kernel(::MyFlux) = my_flux_kernel
```

where `my_flux_kernel(i, j, grid, Tu, clock, fields, parameters)` returns a
heat flux with `parameters.flux::MyFlux`.

Built-in dispatches:

- `ConductiveFlux` → `slab_internal_heat_flux` (single-layer Fourier)
- `IceSnowConductiveFlux` → `ice_snow_conductive_flux` (resistors in series)
- `Function` → returned directly (the function is its own kernel)
"""
@inline flux_kernel(::ConductiveFlux) = slab_internal_heat_flux
@inline flux_kernel(::IceSnowConductiveFlux) = ice_snow_conductive_flux
@inline flux_kernel(f::Function) = f

#####
##### Checkpointing
#####

function prognostic_state(therm::SlabThermodynamics)
    return (top_surface_temperature = prognostic_state(therm.top_surface_temperature),)
end

function restore_prognostic_state!(therm::SlabThermodynamics, state)
    restore_prognostic_state!(therm.top_surface_temperature, state.top_surface_temperature)
    return therm
end
