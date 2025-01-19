using ClimaSeaIce
using ClimaOcean
using ClimaOcean.ECCO
using Adapt
using Oceananigans
using Oceananigans.Operators: intrinsic_vector
using Oceananigans.Grids: node
using Oceananigans.Fields: interpolate, instantiated_location, flattened_unique_values
using Oceananigans.Models: extract_field_time_series
using Oceananigans.Operators

using Oceananigans.OutputReaders: update_field_time_series!, extract_field_time_series

# We extend the τx and τy methods to compute the time-dependent stress
import ClimaSeaIce.SeaIceMomentumEquations: explicit_τx, 
                                            explicit_τy,
                                            implicit_τx_coefficient,
                                            implicit_τy_coefficient

import Oceananigans.Fields: fractional_indices

@inline fractional_indices(at_node, grid::ImmersedBoundaryGrid, ℓx, ℓy, ℓz) = fractional_indices(at_node, grid.underlying_grid, ℓx, ℓy, ℓz) 

#####
##### Defining the types
#####

struct PrescribedOceanStress{FT}
    ρₒ :: FT
    Cᴰ :: FT
end

Adapt.adapt_structure(to, τ::PrescribedOceanStress) = 
    PrescribedOceanStress(Adapt.adapt(to, τ.ρₒ),
                          Adapt.adapt(to, τ.Cᴰ))

@inline explicit_τx(i, j, k, grid, ::PrescribedOceanStress, clock, fields) = zero(grid)
@inline explicit_τy(i, j, k, grid, ::PrescribedOceanStress, clock, fields) = zero(grid)

@inline function implicit_τx_coefficient(i, j, k, grid, τ::PrescribedOceanStress, clock, fields) 
    uᵢ = @inbounds fields.u[i, j, k]
    vᵢ = ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    
    return τ.ρₒ * τ.Cᴰ * sqrt(uᵢ^2 + vᵢ^2)
end

@inline function implicit_τy_coefficient(i, j, k, grid, τ::PrescribedOceanStress, clock, fields) 
    uᵢ = ℑxyᶠᶜᵃ(i, j, k, grid, fields.u)
    vᵢ = @inbounds fields.v[i, j, k]
    
    return τ.ρₒ * τ.Cᴰ * sqrt(uᵢ^2 + vᵢ^2)
end

@inline function interpolate_to_sea_ice_grid(i, j, k, native_grid, fts, grid, time)
    times = fts.times
    data = fts.data
    time_indexing = fts.time_indexing
    backend = fts.backend
    loc = instantiated_location(fts)
    X = node(i, j, k, grid, loc...)

    # Interpolate field time series data onto the current node and time
    return interpolate(X, time, data, loc, native_grid, times, backend, time_indexing)
end    

#####
##### Atmosphere external stress
#####

struct PrescribedAtmosphereStress{U, V, FT, G}
    uₐ :: U
    vₐ :: V
    ρₐ :: FT
    Cᴰ :: FT
    grid :: G
end

Adapt.adapt_structure(to, τ::PrescribedAtmosphereStress) = 
     PrescribedAtmosphereStress(Adapt.adapt(to, τ.uₐ), 
                                Adapt.adapt(to, τ.vₐ), 
                                Adapt.adapt(to, τ.ρₐ),
                                Adapt.adapt(to, τ.Cᴰ),
                                Adapt.adapt(to, τ.grid))

@inline function explicit_τx(i, j, k, grid, τ::PrescribedAtmosphereStress, clock, fields) 
    time = Time(clock.time)
    native_grid = τ.grid

    to_loc   = (Face(),   Center(), nothing)
    from_loc = (Center(), Center(), nothing)

    uₑ, vₑ = interpolate_atmos_velocities(i, j, k, native_grid, to_loc, from_loc, τ.uₐ, τ.vₐ, grid, time)
    
    return τ.ρₐ * τ.Cᴰ * sqrt(uₑ^2 + vₑ^2) * uₑ
end

@inline function explicit_τy(i, j, k, grid, τ::PrescribedAtmosphereStress, clock, fields) 
    time = Time(clock.time)
    native_grid = τ.grid

    to_loc   = (Center(), Face(),   nothing)
    from_loc = (Center(), Center(), nothing)

    uₑ, vₑ = interpolate_atmos_velocities(i, j, k, native_grid, to_loc, from_loc, τ.uₐ, τ.vₐ, grid, time)

    return τ.ρₐ * τ.Cᴰ * sqrt(uₑ^2 + vₑ^2) * vₑ
end

@inline function interpolate_atmos_velocities(i, j, k, native_grid, to_loc, from_loc, u, v, grid, time)
    times = u.times
    u_data = u.data
    v_data = v.data
    time_indexing = u.time_indexing
    backend = u.backend
    X = node(i, j, k, grid, to_loc...)

    uₑ = interpolate(X, time, u_data, from_loc, native_grid, times, backend, time_indexing)
    vₑ = interpolate(X, time, v_data, from_loc, native_grid, times, backend, time_indexing)

    return intrinsic_vector(i, j, k, grid, uₑ, vₑ)
end    