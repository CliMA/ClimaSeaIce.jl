using Oceananigans.Grids: static_column_depthᶜᶜᵃ

"""
    LandfastBasalStress{FT}

Stress exerted by the sea floor on grounded sea-ice keels, following
[Lemieux et al. (2015)](@cite Lemieux2015).

Where the ice is thick enough that its keels reach the sea floor, the bed arrests the ice, producing
the landfast belt observed along shallow shelves and through narrow channels. The stress is

```math
τᵇ = k₂ \\max(0, h - hᶜ) \\exp[-C (1 - ℵ)] \\frac{𝐮}{|𝐮| + u₀}, \\qquad hᶜ = \\frac{H ℵ}{k₁}
```

where ``H`` is the still-water column depth, ``h`` the ice thickness, ``ℵ`` the ice concentration, and
``u₀`` a small velocity that keeps the stress finite as ``𝐮 → 0``. Note that ``τᵇ`` saturates at
``k₂ (h - hᶜ)`` and vanishes with ``𝐮``: the bed can arrest ice, never accelerate it. The stress acts
only where ``H`` is smaller than `maximum_water_depth`.

Because ``τᵇ`` is stiff at small velocity it is treated implicitly by the momentum solver, alongside
the ice-ocean drag. It is carried by the sea floor rather than by the water column, so it is held
separately from `external_momentum_stresses` and never appears in the stress passed to an ocean.

Pass it to [`SeaIceMomentumEquation`](@ref) as `basal_stress`.
"""
struct LandfastBasalStress{FT}
    critical_thickness_parameter :: FT
    stress_parameter :: FT
    concentration_hardening :: FT
    minimum_speed :: FT
    maximum_water_depth :: FT
end

"""
    LandfastBasalStress(FT = Oceananigans.defaults.FloatType;
                        critical_thickness_parameter = 8,
                        stress_parameter = 15,
                        concentration_hardening = 20,
                        minimum_speed = 5e-5,
                        maximum_water_depth = 30)

Construct a `LandfastBasalStress`. The still-water column depth ``H`` is read from the grid, so the
same object may be passed to any grid.

Keyword Arguments
=================

- `critical_thickness_parameter`: ``k₁``, setting the keel thickness that grounds in a column of
                                  depth ``H``. Default: `8`.
- `stress_parameter`: ``k₂`` (N m⁻³). Default: `15`.
- `concentration_hardening`: ``C``, suppressing the stress in unconsolidated ice. Default: `20`.
- `minimum_speed`: ``u₀`` (m s⁻¹), regularizing the stress at vanishing velocity. Default: `5e-5`.
- `maximum_water_depth`: depth (m) beyond which no grounding is possible. Default: `30`.
"""
function LandfastBasalStress(FT::DataType = Oceananigans.defaults.FloatType;
                             critical_thickness_parameter = 8,
                             stress_parameter = 15,
                             concentration_hardening = 20,
                             minimum_speed = 5e-5,
                             maximum_water_depth = 30)

    return LandfastBasalStress(convert(FT, critical_thickness_parameter),
                               convert(FT, stress_parameter),
                               convert(FT, concentration_hardening),
                               convert(FT, minimum_speed),
                               convert(FT, maximum_water_depth))
end

Base.summary(::LandfastBasalStress{FT}) where FT = "LandfastBasalStress{$FT}"

function Base.show(io::IO, b::LandfastBasalStress)
    print(io, summary(b), '\n')
    print(io, "├── critical_thickness_parameter: ", b.critical_thickness_parameter, '\n')
    print(io, "├── stress_parameter: ", b.stress_parameter, '\n')
    print(io, "├── concentration_hardening: ", b.concentration_hardening, '\n')
    print(io, "├── minimum_speed: ", b.minimum_speed, '\n')
    print(io, "└── maximum_water_depth: ", b.maximum_water_depth)
end

Adapt.adapt_structure(to, b::LandfastBasalStress) =
    LandfastBasalStress(b.critical_thickness_parameter,
                        b.stress_parameter,
                        b.concentration_hardening,
                        b.minimum_speed,
                        b.maximum_water_depth)

materialize_basal_stress(stress, grid) = stress

# Grounded keel thickness in excess of what the column can accommodate, hardened by concentration.
@inline function basal_stress_magnitude(i, j, k, grid, b, fields)
    h = @inbounds fields.h[i, j, 1]
    ℵ = @inbounds fields.ℵ[i, j, 1]
    H = static_column_depthᶜᶜᵃ(i, j, grid)
    δh = max(0, h - H * ℵ / b.critical_thickness_parameter)
    kᵇ = b.stress_parameter * δh * exp(- b.concentration_hardening * (1 - ℵ))
    return ifelse(H < b.maximum_water_depth, kᵇ, zero(grid))
end

"""
    basal_τx_coefficient(i, j, k, grid, basal_stress, fields)
    basal_τy_coefficient(i, j, k, grid, basal_stress, fields)

Return the coefficient ``τᵇ / u`` of the basal stress, for implicit treatment by the momentum solver.
"""
@inline function basal_τx_coefficient(i, j, k, grid, b::LandfastBasalStress, fields)
    kᵇ = ℑxᶠᵃᵃ(i, j, 1, grid, basal_stress_magnitude, b, fields)
    u  = @inbounds fields.u[i, j, k]
    v  = ℑxyᶠᶜᵃ(i, j, k, grid, fields.v)
    return kᵇ / (sqrt(u^2 + v^2) + b.minimum_speed)
end

@inline function basal_τy_coefficient(i, j, k, grid, b::LandfastBasalStress, fields)
    kᵇ = ℑyᵃᶠᵃ(i, j, 1, grid, basal_stress_magnitude, b, fields)
    u  = ℑxyᶜᶠᵃ(i, j, k, grid, fields.u)
    v  = @inbounds fields.v[i, j, k]
    return kᵇ / (sqrt(u^2 + v^2) + b.minimum_speed)
end

@inline basal_τx_coefficient(i, j, k, grid, ::Nothing, fields) = zero(grid)
@inline basal_τy_coefficient(i, j, k, grid, ::Nothing, fields) = zero(grid)
