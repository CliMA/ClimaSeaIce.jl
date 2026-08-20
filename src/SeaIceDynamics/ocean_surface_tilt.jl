struct OceanSurfaceTilt{E, S, FT}
    η  :: E   # ocean surface height read by the kernel (an extended copy after materialization)
    η₀ :: S   # source ocean surface height (e.g. a live ocean free surface); copied into η each time step
    g  :: FT
end

"""
    OceanSurfaceTilt(FT = Oceananigans.defaults.FloatType;
                     η = ZeroField(FT),
                     gravitational_acceleration = Oceananigans.defaults.gravitational_acceleration)

A structure representing the acceleration felt by the sea ice because of the tilt of the ocean
surface, computed as

```math
\\begin{align*}
aᵤ & = - g \\frac{∂η}{∂x} \\\\
aᵥ & = - g \\frac{∂η}{∂y}
\\end{align*}
```

where ``η`` is the ocean surface height and ``g`` is the gravitational acceleration.

Arguments
==========
- `FT`: The floating point type of the tilt parameters (optional, default: `Oceananigans.defaults.FloatType`).

Keyword Arguments
==================
- `η`: The ocean surface height, for example the free surface displacement of an ocean model.
- `gravitational_acceleration`: The gravitational acceleration.
"""
function OceanSurfaceTilt(FT = Oceananigans.defaults.FloatType;
                          η = ZeroField(FT),
                          gravitational_acceleration = Oceananigans.defaults.gravitational_acceleration)

    return OceanSurfaceTilt(η, η, convert(FT, gravitational_acceleration))
end

# Drop the source surface height on the device.
Adapt.adapt_structure(to, tilt::OceanSurfaceTilt) =
    OceanSurfaceTilt(Adapt.adapt(to, tilt.η), nothing, tilt.g)

Base.summary(tilt::OceanSurfaceTilt) = string("OceanSurfaceTilt with g = ", tilt.g)

function Base.show(io::IO, tilt::OceanSurfaceTilt)
    print(io, "OceanSurfaceTilt", '\n')
    print(io, "├── η: ", summary(tilt.η), '\n')
    print(io, "└── g: ", tilt.g)
end

materialize_surface_tilt(tilt, grid) = tilt

function materialize_surface_tilt(tilt::OceanSurfaceTilt, grid)
    η = extended_external_variable(tilt.η₀, grid)
    return OceanSurfaceTilt(η, tilt.η₀, tilt.g)
end

# Fill the ocean surface height halos once per time step, before substepping (coupler owns the interior).
update_surface_tilt!(tilt, grid) = nothing

function update_surface_tilt!(tilt::OceanSurfaceTilt, grid)
    tilt.η === tilt.η₀ || refresh_and_fill_external_field!(tilt.η, tilt.η₀)
    return nothing
end

@inline x_surface_tilt(i, j, k, grid, tilt) = zero(grid)
@inline y_surface_tilt(i, j, k, grid, tilt) = zero(grid)

@inline x_surface_tilt(i, j, k, grid, tilt::OceanSurfaceTilt) = tilt.g * ∂xᶠᶜᶜ(i, j, k, grid, tilt.η)
@inline y_surface_tilt(i, j, k, grid, tilt::OceanSurfaceTilt) = tilt.g * ∂yᶜᶠᶜ(i, j, k, grid, tilt.η)
