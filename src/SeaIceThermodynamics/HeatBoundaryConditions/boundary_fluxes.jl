using Oceananigans.Utils: prettysummary

#####
##### Fluxes
#####

# Flux extractor function
# getflux(flux, i, j, grid, T_top, clock, model_fields)

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
        return new{P, T, F}(func, parameters)
    end
end


function Adapt.adapt_structure(to, ff::FluxFunction{P, T}) where {P, T}
    return FluxFunction{T}(adapt(to, ff.func),
                           adapt(to, ff.parameters))
end

Base.summary(flux::FluxFunction{<:Nothing}) = string("FluxFunction of ", prettysummary(flux.func, false))

Base.summary(flux::FluxFunction) = string("FluxFunction of ",
                                          prettysummary(flux.func, false),
                                          " with parameters ",
                                          prettysummary(flux.parameters))

Base.show(io::IO, flux::FluxFunction) = print(io, summary(flux))

"""
    FluxFunction(func; parameters=nothing, top_temperature_dependent=false)

Return `FluxFunction` representing a flux across an air-ice, air-snow, or ice-water interface.
The flux is computed by `func` with the signature

```julia
flux = func(i, j, grid, clock, top_temperature, model_fields)
```

if `isnothing(parameters)`, or

```julia
flux = func(i, j, grid, clock, top_temperature, model_fields, parameters)
```

if `!isnothing(parameters)`. If `func` is `top_temperature_dependent`, then it will be recomputed
during a diagnostic solve for the top temperature.
"""
function FluxFunction(func; parameters=nothing, top_temperature_dependent=false)
    T = top_temperature_dependent ? SurfaceTemperatureDependent() : Nothing
    return FluxFunction{T}(func, parameters)
end

@inline getflux(flux::FluxFunction{<:Nothing}, args...) = flux.func(args...)
@inline getflux(flux::FluxFunction, args...) = flux.func(args..., flux.parameters)

#####
##### Implementations
#####


struct RadiativeEmission{FT}
    emissivity :: FT
    stefan_boltzmann_constant :: FT
    reference_temperature :: FT
end

"""
    RadiativeEmission(FT=Float64; kw...)

Returns a flux representing radiative emission from a surface.
"""
function RadiativeEmission(FT=Float64;
                           emissivity = 1,
                           stefan_boltzmann_constant = 5.67e-8,
                           reference_temperature = 273.15)

    return RadiativeEmission(convert(FT, emissivity),
                             convert(FT, stefan_boltzmann_constant), 
                             convert(FT, reference_temperature))
end

@inline function getflux(emission::RadiativeEmission, i, j, grid, T, clock, fields)
    ϵ = emission.emissivity
    σ = emission.stefan_boltzmann_constant
    Tᵣ = emission.reference_temperature
    return ϵ * σ * (T + Tᵣ)^4 
end

function Base.summary(flux::RadiativeEmission)
    σ = flux.stefan_boltzmann_constant
    ϵ = flux.emissivity
    Tᵣ = flux.reference_temperature
    return string("RadiativeEmission(",
                  "emissivity = ", prettysummary(ϵ), ", ",
                  "stefan_boltzmann_constant = ", prettysummary(σ), ", ",
                  "reference_temperature = ", prettysummary(Tᵣ), ")")
end

Base.show(io::IO, flux::RadiativeEmission) = print(io, summary(flux))

#####
##### show
#####

flux_summary(flux, padchar=" ") = prettysummary(flux)

function flux_summary(fluxes::Tuple, padchar=" ")
    Nfluxes = length(fluxes)
    if Nfluxes == 1
        return string("Tuple with 1 flux:", "\n",
                      "$padchar   └── ", summary(fluxes[1]))
    else
        return string("Tuple with $Nfluxes fluxes:", "\n",
         Tuple(string("$padchar   ├── ", summary(c), "\n") for c in fluxes[1:end-1])...,
                      "$padchar   └── ", summary(fluxes[end]))
    end
end
