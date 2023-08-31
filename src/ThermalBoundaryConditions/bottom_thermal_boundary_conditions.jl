using ClimaSeaIce: melting_temperature

#####
##### Bottom thermal boundary conditions
#####

struct IceWaterThermalEquilibrium{S}
    salinity :: S
end

"""
    IceWaterThermalEquilibrium(; external_fluxes::Tuple=tuple(), salinity=0)

Represents an ice-water interface in thermal equilibrium, such that
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

@inline bottom_temperature(i, j, grid, bottom_bc::PrescribedTemperature, liquidus) = bs.temperature[i, j, 1]

@inline function bottom_temperature(i, j, grid, bottom_bc::IceWaterThermalEquilibrium, liquidus)
    Sₒ = @inbounds bottom_bc.salinity[i, j, 1]
    return melting_temperature(liquidus, Sₒ)
end

@inline function bottom_flux_imbalance(i, j, grid, internal_fluxes, external_fluxes, bottom_thermal_bc, clock, bottom_temperature, model_fields)

    #          
    #   ice        ↑   Qi ≡ internal_fluxes. Example: Qi = - k ∂z T
    #            |⎴⎴⎴|
    # ----------------------- ↕ hᵇ → ℒ ∂t hᵇ = δQ _given_ T = Tₘ
    #            |⎵⎵⎵|
    #   water      ↑   Qx ≡ external_fluxes
    #        

    Qi = getflux(internal_fluxes, i, j, grid, clock, bottom_temperature, model_fields)
    Qx = getflux(external_fluxes, i, j, grid, clock, bottom_temperature, model_fields)

    return Qi - Qx
end


