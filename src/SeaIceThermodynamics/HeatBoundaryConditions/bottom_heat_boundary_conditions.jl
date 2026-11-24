using ...SeaIceThermodynamics: melting_temperature

#####
##### Bottom heat boundary conditions
#####

struct IceWaterThermalEquilibrium{S}
    salinity :: S
end

Adapt.adapt_structure(to, iwte::IceWaterThermalEquilibrium) = IceWaterThermalEquilibrium(adapt(to, iwte.salinity))

"""
    IceWaterThermalEquilibrium(; salinity = 0)

Represents an ice-water interface in heat equilibrium, such that the bottom temperature ``T_b`` is equal to the
melting temperature at the base of the ice,

```math
T_b = Tₘ(S, z_b) ,
```

where ``S`` is the `salinity` at the ice-water boundary and ``z_b`` is the height of the base of floating ice
relative to the sea surface (see `ice_base_height`), set by the ice and snow load and by the `liquid_density`
of the `PhaseTransitions`.

Both freezing and melting may occur at an ice-water boundary.
"""
IceWaterThermalEquilibrium(; salinity = 0) = IceWaterThermalEquilibrium(salinity)

"""
    ice_base_height(hi, hs, ρi, ρs, ρw)

Return the height of the base of floating ice of thickness `hi` and density `ρi`, covered by snow of thickness
`hs` and density `ρs`, relative to the sea surface,

```math
z_b = - \\frac{ρ_i h_i + ρ_s h_s}{ρ_w} ,
```

where ``ρ_w`` is the density of the water. The water pressure at ``z_b`` is the weight of the overlying ice and snow,
independently from the elevation of the free surface.
"""
@inline ice_base_height(hi, hs, ρi, ρs, ρw) = - (ρi * hi + ρs * hs) / ρw

@inline function ice_base_height(i, j, fields, ρw)
    @inbounds begin
        hi = fields.h[i, j, 1]
        ρi = fields.ρi[i, j, 1]
        hs = haskey(fields, :hs) ? fields.hs[i, j, 1] : zero(hi)
        ρs = haskey(fields, :ρs) ? fields.ρs[i, j, 1] : zero(ρi)
    end
    return ice_base_height(hi, hs, ρi, ρs, ρw)
end

@inline bottom_temperature(i, j, grid, bc::PrescribedTemperature, args...) = @inbounds bc.temperature[i, j]
@inline bottom_temperature(i, j, grid, bc::PrescribedTemperature{<:Number}, args...) = bc.temperature

@inline function bottom_temperature(i, j, grid, bc::IceWaterThermalEquilibrium, phase_transitions, fields)
    Sₒ = get_tracer(i, j, 1, grid, bc.salinity)
    zᵇ = ice_base_height(i, j, fields, phase_transitions.liquid_density)
    return melting_temperature(phase_transitions.liquidus, Sₒ, zᵇ)
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
