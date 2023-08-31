#####
##### Fluxes
#####

# Flux extractor function
# getflux(flux, i, j, grid, clock, T_surface, model_fields)

@inline getflux(flux::Number, i, j, grid, args...) = flux
@inline getflux(flux::AbstractArray{<:Any, 2}, i, j, grid, args...) = @inbounds flux[i, j]

# Make sure to check the size on flux construction
@inline getflux(flux::AbstractArray{<:Any, 3}, i, j, grid, args...) = @inbounds flux[i, j, 1]

# Tuples
@inline getflux(fluxes::Tuple{}, i, j, grid, args...) = zero(grid)
@inline getflux(fluxes::Tuple{<:Any}, args...) = getflux(fluxes[1], args...)
@inline getflux(fluxes::Tuple{<:Any, <:Any}, args...) = getflux(fluxes[1], args...) + getflux(fluxes[2], args...)
@inline getflux(fluxes::Tuple, args...) = getflux(fluxes[1], args...) + _getflux(fluxes[2:end], args...)

@inline _getflux(fluxes::Tuple{}, i, j, grid, args...) = zero(grid)
@inline _getflux(fluxes::Tuple{<:Any}, grid, args...) = getflux(fluxes[1], args...)
@inline _getflux(fluxes::Tuple, args...) = getflux(fluxes[1], args...) + getflux(fluxes[2:end], args...)

struct SurfaceTemperatureDependent end

struct FluxFunction{P, T, F}
    func :: F
    parameters :: P
    function FluxFunction{SurfaceTemperatureDependence}(func, parameters=nothing) where SurfaceTemperatureDependence
        T = SurfaceTemperatureDependence
        P = typeof(parameters)
        F = typeof(func)
        new{P, T, F}(func, parameters)
    end
end

"""
    FluxFunction(func; parameters=nothing, surface_temperature_dependent=false)

Return `FluxFunction` representing a flux across an air-ice, air-snow, or ice-water interface.
The flux is computed by `func` with the signature

```julia
flux = func(i, j, grid, clock, surface_temperature, model_fields)
```

if `isnothing(parameters)`, or

```julia
flux = func(i, j, grid, clock, surface_temperature, model_fields, parameters)
```

if `!isnothing(parameters)`. If `func` is `surface_temperature_dependent`, then it will be recomputed
during a diagnostic solve for the surface temperature.
"""
function FluxFunction(func; parameters=nothing, surface_temperature_dependent=false)
    T = surface_temperature_dependent ? SurfaceTemperatureDependent() : Nothing
    return FluxFunction{T}(func, parameters)
end

@inline getflux(flux::FluxFunction{<:Nothing}, args...) = flux.func(args...)
@inline getflux(flux::FluxFunction, args...) = flux.func(args..., flux.parameters)

#####
##### Implementations
#####

@inline _emission(i, j, grid, clock, T, fields, p) = p.ϵ * p.σ * (T + p.Tᵣ)^4

function emission(FT = Float64;
                  emissivity = 1,
                  stefan_boltzmann_constant = 5.67e-8,
                  reference_temperature = 273.15)

    parameters = (ϵ = convert(FT, emissivity),
                  σ = convert(FT, stefan_boltzmann_constant),
                  Tᵣ = convert(FT, reference_temperature))

    return FluxFunction(_emission; parameters, surface_temperature_dependent=true)
end

