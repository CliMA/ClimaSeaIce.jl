using ClimaSeaIce.SeaIceThermodynamics: melting_temperature

#####
##### Bottom heat boundary conditions
#####

struct IceWaterThermalEquilibrium{S}
    salinity :: S
end

Adapt.adapt_structure(to, iwte::IceWaterThermalEquilibrium) =
    IceWaterThermalEquilibrium(adapt(to, iwte.salinity))

"""
    IceWaterThermalEquilibrium(; external_fluxes::Tuple=tuple(), salinity=0)

Represents an ice-water interface in heat equilibrium, such that
the bottom temperature ``T_b`` is equal to the melting temperature,

```math
T_b = Tₘ(S)
```

where ``S`` is the `salinity` at the ice-water boundary.

Both freezing and melting may occur at an ice-water boundary.
The rate of change of the ice-water interface, ``d h_b / dt``,
is related to the sum of `fluxes` into the ice-water interface,

"""
IceWaterThermalEquilibrium(; salinity=0) = IceWaterThermalEquilibrium(salinity)

@inline bottom_temperature(i, j, grid, bc::PrescribedTemperature, args...) = @inbounds bc.temperature[i, j]
@inline bottom_temperature(i, j, grid, bc::PrescribedTemperature{<:Number}, args...) = bc.temperature

@inline function bottom_temperature(i, j, grid, bc::IceWaterThermalEquilibrium, liquidus)
    Sₒ = @inbounds bc.salinity[i, j, 1]
    return melting_temperature(liquidus, Sₒ)
end

@inline function bottom_flux_imbalance(i, j, grid, bottom_heat_bc, top_temperature,
                                       internal_fluxes, external_fluxes, clock, model_fields)

    #
    #   ice        ↑   Qi ≡ internal_fluxes. Example: Qi = - k ∂z T
    #            |⎴⎴⎴|
    # ----------------------- ↕ hᵇ → ℒ ∂t hᵇ = δQ _given_ T = Tₘ
    #            |⎵⎵⎵|
    #   water      ↑   Qx ≡ external_fluxes
    #

    Qi = getflux(internal_fluxes, i, j, grid, top_temperature, clock, model_fields)
    Qx = getflux(external_fluxes, i, j, grid, top_temperature, clock, model_fields)

    return Qi - Qx
end
