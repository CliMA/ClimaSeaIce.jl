# Couple the resolved ice column to a finite slab of snow. Column fluxes are
# upward-positive per ice-covered area; frazil is per grid-cell area. Snow has
# no heat capacity; its latent reservoir is spent before heat reaches the ice.
const snow_conduction_thickness_cutoff = 1e-8 # m; thinner films are thermally indistinguishable from bare ice

struct SnowColumnTopFlux{F, H, I, K, S, M}
    atmosphere :: F
    snow_thickness :: H
    ice_thickness :: I
    snow_conductivity :: K
    itd_shape :: S
    snow_melt_potential :: M
end

Adapt.@adapt_structure SnowColumnTopFlux

@inline function snow_column_conductance(i, j, grid, f::SnowColumnTopFlux)
    @inbounds begin
        h = f.ice_thickness[i, j, 1]
        hs = f.snow_thickness[i, j, 1]
    end
    return f.snow_conductivity * itd_factor(f.itd_shape, h) / max(hs, oftype(hs, snow_conduction_thickness_cutoff))
end

@inline function snow_column_fluxes(i, j, grid, f::SnowColumnTopFlux, Tice, clock, model_fields)
    @inbounds hs = f.snow_thickness[i, j, 1]
    K = snow_column_conductance(i, j, grid, f)
    Ts = column_surface_temperature_balance(i, j, grid, Tice, K, Tice, f.atmosphere, clock, model_fields)
    Qatm = getflux(f.atmosphere, i, j, grid, Ts, clock, model_fields)
    Qmelt = ifelse(Ts >= 0, min(max(K * Tice - Qatm, zero(Tice)),
                              column_parameter(f.snow_melt_potential, i, j)), zero(Tice))
    Qbare = getflux(f.atmosphere, i, j, grid, Tice, clock, model_fields)
    snowy = hs > snow_conduction_thickness_cutoff
    return ifelse(snowy, Ts, Tice), ifelse(snowy, Qatm + Qmelt, Qbare), ifelse(snowy, Qmelt, zero(Tice))
end

@inline function HeatBoundaryConditions.getflux(f::SnowColumnTopFlux, i, j, grid, Tice, clock, model_fields)
    _, Qice, _ = snow_column_fluxes(i, j, grid, f, Tice, clock, model_fields)
    return Qice
end

@inline function column_surface_temperature_balance(i, j, grid, Tₛ⁻, conductance, T_top,
                                                    f::SnowColumnTopFlux, clock, model_fields)
    C = max(conductance, eps(eltype(grid)))
    K = snow_column_conductance(i, j, grid, f)
    # Ice and snow conductances act in series while both surfaces are cold.
    Ts = column_surface_temperature_balance(i, j, grid, Tₛ⁻, C * K / (C + K), T_top,
                                            f.atmosphere, clock, model_fields)
    Ti = (C * T_top + K * Ts) / (C + K)
    Qatm = getflux(f.atmosphere, i, j, grid, zero(T_top), clock, model_fields)
    capacity = column_parameter(f.snow_melt_potential, i, j)
    exhausted = (Ts >= 0) & (K * Ti - Qatm >= capacity)
    Ti = ifelse(exhausted, T_top - (Qatm + capacity) / C, Ti)
    bare = column_surface_temperature_balance(i, j, grid, Tₛ⁻, C, T_top, f.atmosphere, clock, model_fields)
    hs = @inbounds f.snow_thickness[i, j, 1]
    Ti = ifelse(hs > snow_conduction_thickness_cutoff, Ti, bare)
    return ifelse(conductance > eps(eltype(grid)), min(Ti, zero(Ti)), Tₛ⁻)
end

@inline function column_external_flux_temperature_slope(f::SnowColumnTopFlux, i, j, grid, Tice, clock, model_fields)
    hs = @inbounds f.snow_thickness[i, j, 1]
    K = snow_column_conductance(i, j, grid, f)
    Qatm = getflux(f.atmosphere, i, j, grid, zero(Tice), clock, model_fields)
    Qexcess = K * Tice - Qatm
    capacity = column_parameter(f.snow_melt_potential, i, j)
    Ts, _, _ = snow_column_fluxes(i, j, grid, f, Tice, clock, model_fields)
    λ = column_external_flux_temperature_slope(f.atmosphere, i, j, grid, Ts, clock, model_fields)
    melting_slope = ifelse(Qexcess >= capacity, zero(Tice), K)
    snowy_slope = ifelse(Qexcess >= 0, melting_slope, K * λ / (K + λ))
    bare_slope = column_external_flux_temperature_slope(f.atmosphere, i, j, grid, Tice, clock, model_fields)
    return ifelse(hs > snow_conduction_thickness_cutoff, snowy_slope, bare_slope)
end

@kernel function _compute_snow_melt_potential!(potential, grid, snow_thickness, snow_density, latent_heat, Δt)
    i, j = @index(Global, NTuple)
    ρs = column_parameter(snow_density, i, j)
    @inbounds potential[i, j, 1] = ρs * latent_heat * snow_thickness[i, j, 1] / Δt
end

@kernel function _column_snow_volume_update!(ice_thickness, ice_concentration,
                                             snow_thickness, consolidation_thickness,
                                             sea_ice_density, snow_density, snowfall,
                                             phase_transitions, z, fields, auxiliary,
                                             grid, heat_boundary_conditions,
                                             external_heat_fluxes, relation, conductivity,
                                             snow_surface_temperature, mass_fluxes,
                                             clock, model_fields, Δt)
    i, j = @index(Global, NTuple)

    @inbounds begin
        hⁿ = ice_thickness[i, j, 1]
        ℵⁿ = ice_concentration[i, j, 1]
        hsⁿ = snow_thickness[i, j, 1]
        hᶜ = consolidation_thickness[i, j, 1]
        ρi = column_parameter(sea_ice_density, i, j)
        ρs = column_parameter(snow_density, i, j)
        Viⁿ = hⁿ * ℵⁿ
        Vsⁿ = hsⁿ * ℵⁿ

        z.hs⁻[i, j, 1] = z.hsⁿ[i, j, 1]
        z.hb⁻[i, j, 1] = z.hbⁿ[i, j, 1]

        ∂t_V, surface_shift = column_volume_tendency(i, j, grid, hⁿ, ℵⁿ, hᶜ, ρi, phase_transitions,
                                                   fields, auxiliary, heat_boundary_conditions,
                                                   external_heat_fluxes, relation, conductivity,
                                                   clock, model_fields, Δt)

        h⁺, ℵ⁺ = ice_volume_update(ProportionalEvolution(), ∂t_V, hⁿ, ℵⁿ, hᶜ, Δt)

        Tice = auxiliary.surface_temperature[i, j, 1]
        Ts, _, Qmelt = snow_column_fluxes(i, j, grid, external_heat_fluxes.top, Tice, clock, model_fields)
        snow_surface_temperature[i, j, 1] = Ts
        Qatm = getflux(external_heat_fluxes.top.atmosphere, i, j, grid, Ts, clock, model_fields)
        Ksnow = snow_column_conductance(i, j, grid, external_heat_fluxes.top)
        auxiliary.snow_surface_residual[i, j, 1] = ifelse((hsⁿ > snow_conduction_thickness_cutoff) & (Ts < 0),
                                                         abs(Ksnow * (Tice - Ts) - Qatm), zero(Tice))
        snowfall_rate = snow_accumulation(i, j, snowfall, ρs, ℵ⁺, clock)
        Vs⁺ = max(zero(Vsⁿ), Vsⁿ + Δt * (ℵ⁺ * snowfall_rate -
                   ℵⁿ * Qmelt / (ρs * phase_transitions.reference_latent_heat)))
        hs⁺ = ifelse(ℵ⁺ > 0, Vs⁺ / ℵ⁺, zero(hsⁿ))

        h_before_flooding = h⁺
        h⁺, hs⁺ = snow_ice_formation(h⁺, hs⁺, ρi, ρs, phase_transitions.liquid_density)
        surface_shift += h⁺ - h_before_flooding
        hs⁺ = ifelse(ℵ⁺ > 0, hs⁺, zero(hs⁺))

        H = ifelse(h⁺ > 0, h⁺, z.hs⁻[i, j, 1] - z.hb⁻[i, j, 1])
        z.hsⁿ[i, j, 1] = z.hs⁻[i, j, 1] + surface_shift
        z.hbⁿ[i, j, 1] = z.hsⁿ[i, j, 1] - H
        ice_thickness[i, j, 1] = h⁺
        ice_concentration[i, j, 1] = ℵ⁺
        snow_thickness[i, j, 1] = hs⁺

        Psᵃᵇˢ = ρs * snowfall_rate * ℵ⁺
        mass_fluxes.thermodynamics.ice[i, j, 1] = ρi * (h⁺ * ℵ⁺ - Viⁿ) / Δt
        mass_fluxes.thermodynamics.snow[i, j, 1] = ρs * (hs⁺ * ℵ⁺ - Vsⁿ) / Δt - Psᵃᵇˢ
        mass_fluxes.intercepted_snowfall[i, j, 1] = Psᵃᵇˢ
    end
end

function thermodynamic_time_step!(model, th::ColumnEnergyThermodynamics,
                                  snow::SlabThermodynamics, Δt)
    apply_column_thermodynamics_this_stage(model.timestepper, model.clock) || return nothing
    th.heat_boundary_conditions.top isa MeltingConstrainedFluxBalance ||
        throw(ArgumentError("Column ice with slab snow requires MeltingConstrainedFluxBalance at the ice surface"))
    snow.heat_boundary_conditions.top isa MeltingConstrainedFluxBalance ||
        throw(ArgumentError("Column ice with slab snow requires MeltingConstrainedFluxBalance at the snow surface"))
    snow.internal_heat_flux isa ConductiveFlux ||
        throw(ArgumentError("Column ice with slab snow requires conductive slab snow"))
    snow.top_surface_temperature === th.auxiliary.surface_temperature &&
        throw(ArgumentError("Snow and ice surface temperatures must use separate fields"))

    cached = cached_external_heat_fluxes(th, model.external_heat_fluxes, model.clock, fields(model))
    grid = model.grid
    launch!(architecture(grid), grid, :xy, _compute_snow_melt_potential!,
            th.auxiliary.snow_melt_potential, grid, model.snow_thickness, model.snow_density,
            model.phase_transitions.reference_latent_heat, Δt)
    snow_top = SnowColumnTopFlux(cached.top, model.snow_thickness, model.ice_thickness,
                                 snow.internal_heat_flux.conductivity,
                                 snow.internal_heat_flux.itd_shape, th.auxiliary.snow_melt_potential)
    Qe = (top = snow_top, bottom = cached.bottom)
    column_energy_time_step!(th, Qe, model.clock, fields(model), model.ice_consolidation_thickness, Δt)

    launch!(architecture(grid), grid, :xy, _column_snow_volume_update!,
            model.ice_thickness, model.ice_concentration, model.snow_thickness,
            model.ice_consolidation_thickness, model.sea_ice_density, model.snow_density,
            model.snowfall, model.phase_transitions, sea_ice_discretization(grid),
            th.fields, th.auxiliary, grid, th.heat_boundary_conditions, Qe,
            th.relation, thermal_conductivity(th.energy_transport),
            snow.top_surface_temperature, model.mass_fluxes,
            model.clock, fields(model), Δt)
    return nothing
end
