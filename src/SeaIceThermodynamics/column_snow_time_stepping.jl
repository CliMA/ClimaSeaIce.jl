# Couple the resolved ice column to a finite slab of snow. Column fluxes are
# upward-positive per ice-covered area; frazil is per grid-cell area. Snow has
# no heat capacity; its latent reservoir is spent before heat reaches the ice.
const snow_conduction_thickness_cutoff = 1e-8 # m; thinner films are thermally indistinguishable from bare ice

struct SnowColumnTopFlux{F, H, I, K, S, D, L, T, DT}
    atmosphere :: F
    snow_thickness :: H
    ice_thickness :: I
    snow_conductivity :: K
    itd_shape :: S
    snow_density :: D
    latent_heat :: L
    snow_surface_temperature :: T
    Δt :: DT
end

Adapt.adapt_structure(to, f::SnowColumnTopFlux) =
    SnowColumnTopFlux(Adapt.adapt(to, f.atmosphere),
                      Adapt.adapt(to, f.snow_thickness),
                      Adapt.adapt(to, f.ice_thickness),
                      Adapt.adapt(to, f.snow_conductivity),
                      Adapt.adapt(to, f.itd_shape),
                      Adapt.adapt(to, f.snow_density),
                      Adapt.adapt(to, f.latent_heat),
                      Adapt.adapt(to, f.snow_surface_temperature), f.Δt)

@inline function snow_column_conductance(f::SnowColumnTopFlux, i, j, hs)
    @inbounds h = f.ice_thickness[i, j, 1]
    return f.snow_conductivity * itd_factor(f.itd_shape, h) / hs
end

@inline function snow_column_fluxes(f::SnowColumnTopFlux, i, j, grid, Tice, clock, model_fields)
    @inbounds hs = f.snow_thickness[i, j, 1]
    if hs <= snow_conduction_thickness_cutoff
        Q = getflux(f.atmosphere, i, j, grid, Tice, clock, model_fields)
        return Tice, Q, zero(Q)
    end

    K = snow_column_conductance(f, i, j, hs)
    balance(Ts) = K * (Tice - Ts) - getflux(f.atmosphere, i, j, grid, Ts, clock, model_fields)
    tolerance(Ts) = column_flux_balance_tolerance(grid,
                                                   getflux(f.atmosphere, i, j, grid, Ts, clock, model_fields),
                                                   K * (Tice - Ts), K,
                                                   max(abs(Tice), abs(Ts)))

    # A positive residual at the snow melting point means that the unconstrained
    # skin would be above 0 C. Pin it at melting and spend the excess on snow.
    Qatmₘ = getflux(f.atmosphere, i, j, grid, zero(Tice), clock, model_fields)
    Qconductₘ = K * Tice
    if Qconductₘ >= Qatmₘ
        ρs = column_parameter(f.snow_density, i, j)
        capacity = ρs * f.latent_heat * hs / f.Δt
        Qmelt = min(Qconductₘ - Qatmₘ, capacity)
        return zero(Tice), Qatmₘ + Qmelt, Qmelt
    end

    # On the cold branch, solve snow conduction = atmospheric heat loss.
    # Include the ice temperature in the bracket: a fixed -100°C endpoint
    # does not enclose the root if an upstream ice-surface solve fails cold.
    cold = min(oftype(Tice, -100), Tice - one(Tice))
    warm = zero(Tice)
    width = one(Tice)
    cold_residual = balance(cold)
    for _ in 1:32
        isfinite(cold_residual) || return oftype(Tice, NaN), oftype(Tice, NaN), oftype(Tice, NaN)
        cold_residual >= 0 && break
        width *= 2
        cold -= width
        cold_residual = balance(cold)
    end
    if !(isfinite(cold_residual) && cold_residual >= 0)
        bad = oftype(Tice, NaN)
        return bad, bad, bad
    end

    for _ in 1:64
        trial = (cold + warm) / 2
        residual = balance(trial)
        if !isfinite(residual)
            bad = oftype(Tice, NaN)
            return bad, bad, bad
        end
        abs(residual) <= tolerance(trial) && return trial, getflux(f.atmosphere, i, j, grid, trial, clock, model_fields), zero(Tice)
        if residual > 0
            cold = trial
        else
            warm = trial
        end
    end
    Ts = (cold + warm) / 2
    if abs(balance(Ts)) > tolerance(Ts)
        bad = oftype(Tice, NaN)
        return bad, bad, bad
    end
    Qice = getflux(f.atmosphere, i, j, grid, Ts, clock, model_fields)
    return Ts, Qice, zero(Qice)
end

@inline function HeatBoundaryConditions.getflux(f::SnowColumnTopFlux, i, j, grid, Tice, clock, model_fields)
    _, Qice, _ = snow_column_fluxes(f, i, j, grid, Tice, clock, model_fields)
    return Qice
end

# Differentiate the snow/atmosphere balance with respect to the ice surface
# temperature. The resulting conductance enters the implicit ice-column row.
@inline function column_external_flux_temperature_slope(f::SnowColumnTopFlux, i, j, grid, Tice, clock, model_fields)
    hs = @inbounds f.snow_thickness[i, j, 1]
    hs <= snow_conduction_thickness_cutoff &&
        return column_external_flux_temperature_slope(f.atmosphere, i, j, grid, Tice, clock, model_fields)

    K = snow_column_conductance(f, i, j, hs)
    Qatmₘ = getflux(f.atmosphere, i, j, grid, zero(Tice), clock, model_fields)
    Qexcess = K * Tice - Qatmₘ
    if Qexcess >= 0
        ρs = column_parameter(f.snow_density, i, j)
        capacity = ρs * f.latent_heat * hs / f.Δt
        return Qexcess >= capacity ? zero(Tice) : K
    end

    Ts, _, _ = snow_column_fluxes(f, i, j, grid, Tice, clock, model_fields)
    λ = column_external_flux_temperature_slope(f.atmosphere, i, j, grid, Ts, clock, model_fields)
    return K * λ / (K + λ)
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
    Nz = size(grid, 3)

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

        if hⁿ >= hᶜ
            Qᵇ = column_basal_stefan_flux(heat_boundary_conditions.bottom,
                                          external_heat_fluxes.bottom, fields, auxiliary,
                                          grid, relation, clock, model_fields, i, j)
            Qˢ = column_surface_stefan_residual_flux(heat_boundary_conditions.top,
                                                     external_heat_fluxes.top, i, j, Nz,
                                                     grid, auxiliary, fields, relation,
                                                     clock, model_fields, Δt)
            δhᵇ = column_stefan_thickness_change(phase_transitions, ρi, Qᵇ, Δt)
            δhˢ = column_surface_melt_thickness_change(Qˢ, Δt, fields, relation, grid, i, j)
            cell_flux_correction = column_cell_flux_correction(external_heat_fluxes.bottom, ℵⁿ, i, j, grid,
                                                               zero(eltype(grid)), clock, model_fields)
            ∂t_V = ℵⁿ * (δhᵇ + δhˢ) / Δt - cell_flux_correction / (ρi * phase_transitions.reference_latent_heat)
            surface_shift = δhˢ
        else
            Tb = column_bottom_reference_temperature(heat_boundary_conditions.bottom,
                                                     i, j, grid, fields, relation)
            ℰ = ρi * latent_heat(phase_transitions, Tb)
            Qu = column_unconsolidated_top_flux(heat_boundary_conditions.top,
                                                external_heat_fluxes.top, conductivity,
                                                Tb, hⁿ, i, j, grid, fields, relation,
                                                clock, model_fields)
            Qb = getflux(external_heat_fluxes.bottom, i, j, grid, Tb, clock, model_fields)
            cell_flux_correction = column_cell_flux_correction(external_heat_fluxes.bottom, ℵⁿ, i, j, grid, Tb, clock, model_fields)
            ∂t_V = (ℵⁿ * (Qu - Qb) - cell_flux_correction) / ℰ
            surface_shift = zero(eltype(grid))
        end

        h⁺, ℵ⁺ = ice_volume_update(ProportionalEvolution(), ∂t_V, hⁿ, ℵⁿ, hᶜ, Δt)

        Tice = auxiliary.surface_temperature[i, j, 1]
        Ts, _, Qmelt = snow_column_fluxes(external_heat_fluxes.top, i, j,
                                          grid, Tice, clock, model_fields)
        snow_surface_temperature[i, j, 1] = Ts
        if hsⁿ > snow_conduction_thickness_cutoff && Ts < 0
            Qatm = getflux(external_heat_fluxes.top.atmosphere, i, j, grid, Ts, clock, model_fields)
            Ksnow = snow_column_conductance(external_heat_fluxes.top, i, j, hsⁿ)
            auxiliary.snow_surface_residual[i, j, 1] = abs(Ksnow * (Tice - Ts) - Qatm)
        else
            auxiliary.snow_surface_residual[i, j, 1] = zero(Tice)
        end
        snowfall_rate = snow_accumulation(i, j, snowfall, ρs, ℵ⁺, clock)
        hs_conservative = ifelse(ℵ⁺ > 0, Vsⁿ / ℵ⁺, zero(hsⁿ))
        hs⁺ = max(zero(hsⁿ), hs_conservative +
                  Δt * (snowfall_rate - Qmelt / (ρs * phase_transitions.reference_latent_heat)))

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
    snow_top = SnowColumnTopFlux(cached.top, model.snow_thickness, model.ice_thickness,
                                 snow.internal_heat_flux.conductivity,
                                 snow.internal_heat_flux.itd_shape,
                                 model.snow_density,
                                 model.phase_transitions.reference_latent_heat,
                                 snow.top_surface_temperature, Δt)
    Qe = (top = snow_top, bottom = cached.bottom)
    column_energy_time_step!(th, Qe, model.clock, fields(model), model.ice_consolidation_thickness, Δt)

    grid = model.grid
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
