"""
    ColumnBasalHeatFlux(interface, interface_cell, frazil)

Upward-positive basal fluxes for a coupled ice column. `interface` is per ice area;
`interface_cell` and `frazil` are per grid-cell area.
"""
struct ColumnBasalHeatFlux{I, C, F}
    interface :: I
    interface_cell :: C
    frazil :: F
end

Adapt.@adapt_structure ColumnBasalHeatFlux

@inline HeatBoundaryConditions.getflux(flux::ColumnBasalHeatFlux, i, j, grid, T, clock, model_fields) =
    getflux(flux.interface, i, j, grid, T, clock, model_fields)

@inline column_cell_flux_correction(i, j, grid, flux, ice_fraction, T, clock, model_fields) = zero(eltype(grid))
@inline function column_cell_flux_correction(i, j, grid, flux::ColumnBasalHeatFlux, ice_fraction, T, clock, model_fields)
    frazil_cell = getflux(flux.frazil, i, j, grid, T, clock, model_fields)
    interface_cell = getflux(flux.interface_cell, i, j, grid, T, clock, model_fields)
    interface_ice = getflux(flux.interface, i, j, grid, T, clock, model_fields)
    # ice_fraction * interface_ice + correction = interface_cell + frazil_cell.
    return frazil_cell + interface_cell - ice_fraction * interface_ice
end

"""
    column_stefan_thickness_change(phase_transitions, sea_ice_density, residual_energy_flux, Δt)

Return the Stefan thickness change `Δt Q / (ρᵢ ℒ)` implied by a residual interface energy flux `Q`, positive for growth.
"""
@inline column_stefan_thickness_change(𝒫, ρᵢ, Q, Δt) = Δt * Q / (ρᵢ * 𝒫.reference_latent_heat)

# Icepack's surface melt consumes the enthalpy of each layer that disappears,
# including its remaining sensible and brine-pocket latent heat. In this
# shifted-energy convention, -q_ice = E_melt - E + ρᵢ c_w μ S.
@inline function column_melt_energy_density(relation, E, S)
    meltwater_sensible_heat = liquid_volumetric_heat_capacity(relation) * liquidus_slope(relation) * S
    return max(zero(E), complete_melt_energy(relation, S) - E + meltwater_sensible_heat)
end

@inline function column_surface_melt_thickness_change(i, j, grid, Qˢ, Δt, fields, relation)
    remaining_energy = max(zero(Qˢ), -Qˢ * Δt)
    δhˢ = zero(Qˢ)
    Nz = size(grid, 3)
    for k in Nz:-1:1
        @inbounds begin
            E = fields.internal_energy[i, j, k]
            S = fields.bulk_salinity[i, j, k]
        end
        Δz = Δzᶜᶜᶜ(i, j, k, grid)
        energy_density = column_melt_energy_density(relation, E, S)
        melt = ifelse(remaining_energy > 0, min(Δz, remaining_energy / max(energy_density, eps(energy_density))), zero(Δz))
        δhˢ -= melt
        remaining_energy -= melt * energy_density
    end
    return δhˢ
end

# Upward-positive basal Stefan flux. A Dirichlet base grows ice from the conductive minus the ocean heat flux.
@inline column_basal_stefan_flux(i, j, grid, boundary, ext, fields, auxiliary, relation, clock, model_fields) = zero(eltype(grid))

@inline function column_basal_stefan_flux(i, j, grid, boundary::ColumnDirichletBoundary, ext, fields, auxiliary, relation, clock, model_fields)
    conductance = column_boundary_temperature_conductance(i, j, 1, 1, grid, auxiliary)
    T₁ = @inbounds fields.temperature[i, j, 1]
    Tᵇ = column_dirichlet_temperature(boundary, i, j, grid, relation)
    Qᵇ = getflux(ext, i, j, grid, Tᵇ, clock, model_fields)
    return conductance * (Tᵇ - T₁) - Qᵇ
end

# Unconsolidated (single-layer slab) top heat loss; a Dirichlet top adds the single-layer conduction k/h (Tᵇ − Tᵗ).
@inline function column_unconsolidated_top_flux(i, j, grid, bc, ext, conductivity, Tᵇ, h, fields, auxiliary, relation, clock, model_fields)
    Q = getflux(ext, i, j, grid, Tᵇ, clock, model_fields)
    @inbounds surface_flux = auxiliary.surface_fluxes.heat[i, j, 1] - auxiliary.surface_fluxes.ice_melt[i, j, 1]
    return ifelse(bc isa MeltingConstrainedFluxBalance, surface_flux, Q)
end

@inline function column_unconsolidated_top_flux(i, j, grid, bc::ColumnDirichletBoundary, ext, conductivity, Tᵇ, h, fields, auxiliary, relation, clock, model_fields)
    Tᵗ = column_dirichlet_temperature(bc, i, j, grid, relation)
    k  = ice_thermal_conductivity(conductivity, Tᵗ, @inbounds(fields.bulk_salinity[i, j, 1]))
    Qi = ifelse(h ≤ zero(h), zero(h), k / h * (Tᵇ - Tᵗ))
    return Qi + getflux(ext, i, j, grid, Tᵗ, clock, model_fields)
end

@inline function column_volume_tendency(i, j, grid, h, ℵ, hᶜ, ρ, phase_transitions,
                                        fields, auxiliary, bcs, Qe, relation, conductivity,
                                        clock, model_fields, Δt)
    Qᵇ = column_basal_stefan_flux(i, j, grid, bcs.bottom, Qe.bottom, fields, auxiliary, relation, clock, model_fields)
    Qˢ = -auxiliary.surface_fluxes.ice_melt[i, j, 1]
    δhᵇ = column_stefan_thickness_change(phase_transitions, ρ, Qᵇ, Δt)
    δhˢ = column_surface_melt_thickness_change(i, j, grid, Qˢ, Δt, fields, relation)
    correction = column_cell_flux_correction(i, j, grid, Qe.bottom, ℵ, zero(eltype(grid)), clock, model_fields)
    consolidated_tendency = ℵ * (δhᵇ + δhˢ) / Δt - correction / (ρ * phase_transitions.reference_latent_heat)

    Tᵇ = column_bottom_reference_temperature(bcs.bottom, i, j, grid, fields, relation)
    ℰ = ρ * latent_heat(phase_transitions, Tᵇ)
    Qᵘ = column_unconsolidated_top_flux(i, j, grid, bcs.top, Qe.top, conductivity, Tᵇ, h, fields, auxiliary, relation, clock, model_fields)
    Qᵇ = getflux(Qe.bottom, i, j, grid, Tᵇ, clock, model_fields)
    correction = column_cell_flux_correction(i, j, grid, Qe.bottom, ℵ, Tᵇ, clock, model_fields)
    slab_tendency = (ℵ * (Qᵘ - Qᵇ) - correction) / ℰ
    consolidated = column_consolidated(h, hᶜ)
    return ifelse(consolidated, consolidated_tendency, slab_tendency), ifelse(consolidated, δhˢ, zero(δhˢ))
end

@inline column_snow_update!(i, j, grid, ::Nothing, h, ℵ, ℵⁿ, auxiliary, phase_transitions, ρi, snowfall, mass_fluxes, clock, Δt) = (h, zero(h))

@inline function column_snow_update!(i, j, grid, snow, h, ℵ, ℵⁿ, auxiliary, phase_transitions, ρi, snowfall, mass_fluxes, clock, Δt)
    @inbounds begin
        hsⁿ = snow.thickness[i, j, 1]
        ρs = get_tracer(i, j, 1, grid, snow.density)
        Qs = auxiliary.surface_fluxes.snow_melt[i, j, 1]
        Vsⁿ = ℵⁿ * hsⁿ
        Ps = ρs * snow_accumulation(i, j, snowfall, ρs, ℵ, clock) * ℵ
        Vs = max(0, Vsⁿ + Δt * (Ps / ρs - ℵⁿ * Qs / (ρs * phase_transitions.reference_latent_heat)))
        hs = ifelse(ℵ > 0, Vs / ℵ, zero(Vs))
        h⁺, hs⁺ = snow_ice_formation(h, hs, ρi, ρs, phase_transitions.liquid_density)
        hs⁺ = ifelse(ℵ > 0, hs⁺, zero(hs⁺))
        snow.thickness[i, j, 1] = hs⁺
        mass_fluxes.thermodynamics.snow[i, j, 1] = ρs * (hs⁺ * ℵ - Vsⁿ) / Δt - Ps
        mass_fluxes.intercepted_snowfall[i, j, 1] = Ps
    end
    return h⁺, h⁺ - h
end

@kernel function _column_stefan_volume_update!(ice_thickness, ice_concentration, consolidation_thickness,
                                               sea_ice_density, phase_transitions, z, fields, auxiliary, grid,
                                               heat_boundary_conditions, external_heat_fluxes, relation, conductivity,
                                               clock, model_fields, mass_fluxes, snow, snowfall, Δt)
    i, j = @index(Global, NTuple)
    bcs = heat_boundary_conditions
    Qe = external_heat_fluxes

    @inbounds begin
        hⁿ = ice_thickness[i, j, 1]
        ℵⁿ = ice_concentration[i, j, 1]
        hᶜ = consolidation_thickness[i, j, 1]
        ρᵢ = get_tracer(i, j, 1, grid, sea_ice_density)
        Viⁿ = hⁿ * ℵⁿ
        mass_fluxes.thermodynamics.snow[i, j, 1] = 0
        mass_fluxes.intercepted_snowfall[i, j, 1] = 0

        z.hs⁻[i, j, 1] = z.hsⁿ[i, j, 1]
        z.hb⁻[i, j, 1] = z.hbⁿ[i, j, 1]

        ∂t_V, δzˢ = column_volume_tendency(i, j, grid, hⁿ, ℵⁿ, hᶜ, ρᵢ, phase_transitions,
                                          fields, auxiliary, bcs, Qe, relation, conductivity,
                                          clock, model_fields, Δt)

        h⁺, ℵ⁺ = ice_volume_update(ProportionalEvolution(), ∂t_V, hⁿ, ℵⁿ, hᶜ, Δt)
        h⁺, flooding = column_snow_update!(i, j, grid, snow, h⁺, ℵ⁺, ℵⁿ, auxiliary, phase_transitions, ρᵢ, snowfall, mass_fluxes, clock, Δt)
        δzˢ += flooding

        # An ice-free column keeps its previous height so that Δz never vanishes.
        H = ifelse(h⁺ > 0, h⁺, z.hs⁻[i, j, 1] - z.hb⁻[i, j, 1])
        z.hsⁿ[i, j, 1] = z.hs⁻[i, j, 1] + δzˢ
        z.hbⁿ[i, j, 1] = z.hsⁿ[i, j, 1] - H

        ice_thickness[i, j, 1] = h⁺
        ice_concentration[i, j, 1] = ℵ⁺
        mass_fluxes.thermodynamics.ice[i, j, 1] = ρᵢ * (h⁺ * ℵ⁺ - Viⁿ) / Δt
    end
end

function column_stefan_volume_update!(thermodynamics::ColumnEnergyThermodynamics, model, external_heat_fluxes, Δt; snow = nothing)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xy, _column_stefan_volume_update!,
            model.ice_thickness, model.ice_concentration, model.ice_consolidation_thickness,
            model.sea_ice_density, model.phase_transitions, sea_ice_discretization(grid),
            thermodynamics.fields, thermodynamics.auxiliary, grid,
            thermodynamics.heat_boundary_conditions, external_heat_fluxes, thermodynamics.relation,
            thermal_conductivity(thermodynamics.energy_transport), model.clock, fields(model),
            model.mass_fluxes, snow, model.snowfall, Δt)
    return nothing
end
