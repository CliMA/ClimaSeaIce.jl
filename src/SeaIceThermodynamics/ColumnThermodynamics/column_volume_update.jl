@inline column_parameter(value::Number, i, j) = value
@inline column_parameter(value::AbstractField, i, j) = @inbounds value[i, j, 1]

"""
    column_stefan_thickness_change(phase_transitions, sea_ice_density, residual_energy_flux, Δt)

Return the Stefan thickness change `Δt Q / (ρᵢ ℒ)` implied by a residual interface energy flux `Q`, positive for growth.
"""
@inline column_stefan_thickness_change(𝒫, ρᵢ, Q, Δt) = Δt * Q / (ρᵢ * 𝒫.reference_latent_heat)

# Upward-positive basal Stefan flux. A Dirichlet base grows ice from the conductive minus the ocean heat flux.
@inline column_basal_stefan_flux(boundary, ext, fields, auxiliary, grid, relation, clock, model_fields, i, j) = zero(eltype(grid))

@inline function column_basal_stefan_flux(boundary::ColumnDirichletBoundary, ext, fields, auxiliary, grid, relation, clock, model_fields, i, j)
    conductance = column_boundary_temperature_conductance(i, j, 1, 1, grid, auxiliary)
    T₁ = @inbounds fields.temperature[i, j, 1]
    Tᵇ = column_dirichlet_temperature(boundary, i, j, grid, relation)
    Qᵇ = getflux(ext, i, j, grid, Tᵇ, clock, model_fields)
    return conductance * (Tᵇ - T₁) - Qᵇ
end

# Unconsolidated (single-layer slab) top heat loss; a Dirichlet top adds the single-layer conduction k/h (Tᵇ − Tᵗ).
@inline column_unconsolidated_top_flux(bc, ext, conductivity, Tᵇ, h, i, j, grid, fields, relation, clock, model_fields) = getflux(ext, i, j, grid, Tᵇ, clock, model_fields)

@inline function column_unconsolidated_top_flux(bc::ColumnDirichletBoundary, ext, conductivity, Tᵇ, h, i, j, grid, fields, relation, clock, model_fields)
    Tᵗ = column_dirichlet_temperature(bc, i, j, grid, relation)
    k  = ice_thermal_conductivity(conductivity, Tᵗ, @inbounds(fields.bulk_salinity[i, j, 1]))
    Qi = ifelse(h ≤ zero(h), zero(h), k / h * (Tᵇ - Tᵗ))
    return Qi + getflux(ext, i, j, grid, Tᵗ, clock, model_fields)
end

@kernel function _column_stefan_volume_update!(ice_thickness, ice_concentration, consolidation_thickness,
                                               sea_ice_density, phase_transitions, z, fields, auxiliary, grid,
                                               heat_boundary_conditions, external_heat_fluxes, relation, conductivity,
                                               clock, model_fields, Δt)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)
    bcs = heat_boundary_conditions
    Qe = external_heat_fluxes

    @inbounds begin
        hⁿ = ice_thickness[i, j, 1]
        ℵⁿ = ice_concentration[i, j, 1]
        hᶜ = consolidation_thickness[i, j, 1]
        ρᵢ = column_parameter(sea_ice_density, i, j)

        z.hs⁻[i, j, 1] = z.hsⁿ[i, j, 1]
        z.hb⁻[i, j, 1] = z.hbⁿ[i, j, 1]

        if hⁿ ≥ hᶜ
            Qᵇ = column_basal_stefan_flux(bcs.bottom, Qe.bottom, fields, auxiliary, grid, relation, clock, model_fields, i, j)
            Qˢ = column_surface_stefan_residual_flux(bcs.top, Qe.top, i, j, Nz, grid, auxiliary, fields, relation, clock, model_fields, Δt)
            δhᵇ = column_stefan_thickness_change(phase_transitions, ρᵢ, Qᵇ, Δt)
            δhˢ = column_stefan_thickness_change(phase_transitions, ρᵢ, Qˢ, Δt)
            ∂t_V = (δhᵇ + δhˢ) / Δt
            δzˢ = δhˢ
        else
            Tᵇ = column_bottom_reference_temperature(bcs.bottom, i, j, grid, fields, relation)
            ℰ  = ρᵢ * latent_heat(phase_transitions, Tᵇ)
            Qᵘ = column_unconsolidated_top_flux(bcs.top, Qe.top, conductivity, Tᵇ, hⁿ, i, j, grid, fields, relation, clock, model_fields)
            Qᵇ = getflux(Qe.bottom, i, j, grid, Tᵇ, clock, model_fields)
            ∂t_V = (Qᵘ - Qᵇ) / ℰ
            δzˢ = zero(eltype(grid))
        end

        h⁺, ℵ⁺ = ice_volume_update(ProportionalEvolution(), ∂t_V, hⁿ, ℵⁿ, hᶜ, Δt)

        # An ice-free column keeps its previous height so that Δz never vanishes.
        H = ifelse(h⁺ > 0, h⁺, z.hs⁻[i, j, 1] - z.hb⁻[i, j, 1])
        z.hsⁿ[i, j, 1] = z.hs⁻[i, j, 1] + δzˢ
        z.hbⁿ[i, j, 1] = z.hsⁿ[i, j, 1] - H

        ice_thickness[i, j, 1] = h⁺
        ice_concentration[i, j, 1] = ℵ⁺
    end
end

function column_stefan_volume_update!(thermodynamics::ColumnEnergyThermodynamics, model, external_heat_fluxes, Δt)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xy, _column_stefan_volume_update!,
            model.ice_thickness, model.ice_concentration, model.ice_consolidation_thickness,
            model.sea_ice_density, model.phase_transitions, sea_ice_discretization(grid),
            thermodynamics.fields, thermodynamics.auxiliary, grid,
            thermodynamics.heat_boundary_conditions, external_heat_fluxes, thermodynamics.relation,
            thermal_conductivity(thermodynamics.energy_transport), model.clock, fields(model), Δt)
    return nothing
end
