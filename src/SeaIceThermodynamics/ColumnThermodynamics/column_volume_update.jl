@inline column_stefan_parameter(value::Number, i, j, k) = value
@inline column_stefan_parameter(value::AbstractField, i, j, k) = @inbounds value[i, j, k]

"""
    column_stefan_thickness_change(phase_transitions, sea_ice_density, residual_energy_flux, dt)

Return the conservative Stefan thickness change implied by a residual interface
energy flux. `residual_energy_flux` is positive for ice growth.
"""
@inline function column_stefan_thickness_change(phase_transitions,
                                                sea_ice_density,
                                                residual_energy_flux,
                                                Δt)
    volumetric_latent_heat = sea_ice_density * phase_transitions.reference_latent_heat
    return Δt * residual_energy_flux / volumetric_latent_heat
end

@kernel function _column_stefan_thickness_update!(ice_thickness,
                                                  phase_transitions,
                                                  sea_ice_density,
                                                  residual_energy_flux,
                                                  Δt)
    i, j = @index(Global, NTuple)

    @inbounds begin
        ρi = column_stefan_parameter(sea_ice_density, i, j, 1)
        δQ = column_stefan_parameter(residual_energy_flux, i, j, 1)
        ice_thickness[i, j, 1] += column_stefan_thickness_change(phase_transitions, ρi, δQ, Δt)
    end
end

"""
    column_stefan_thickness_update!(ice_thickness, phase_transitions, sea_ice_density, residual_energy_flux, dt)

Update `ice_thickness` conservatively from a residual interface energy flux.
The residual flux is positive for growth and negative for melt.
"""
function column_stefan_thickness_update!(ice_thickness,
                                         phase_transitions,
                                         sea_ice_density,
                                         residual_energy_flux,
                                         Δt)
    grid = ice_thickness.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy,
            _column_stefan_thickness_update!,
            ice_thickness,
            phase_transitions,
            sea_ice_density,
            residual_energy_flux,
            Δt)

    return nothing
end

"""
    column_stefan_thickness_budget(initial_thickness, final_thickness, phase_transitions, sea_ice_density, residual_energy_flux, dt)

Return a named tuple summarizing the conservative Stefan thickness update for
a scalar residual flux.
"""
function column_stefan_thickness_budget(initial_thickness,
                                        final_thickness,
                                        phase_transitions,
                                        sea_ice_density,
                                        residual_energy_flux,
                                        Δt)
    expected_change = column_stefan_thickness_change(phase_transitions,
                                                     sea_ice_density,
                                                     residual_energy_flux,
                                                     Δt)
    actual_change = final_thickness - initial_thickness
    residual = actual_change - expected_change

    return (initial = initial_thickness,
            final = final_thickness,
            actual_change,
            expected_change,
            residual,
            relative_residual = budget_relative_residual(residual,
                                                         actual_change,
                                                         expected_change))
end

function check_column_remap_inputs(source_values, source_faces, target_faces)
    length(source_values) >= 1 ||
        throw(ArgumentError("source_values must contain at least one entry"))

    length(source_faces) == length(source_values) + 1 ||
        throw(ArgumentError("source_faces must have one more entry than source_values"))

    length(target_faces) >= 2 ||
        throw(ArgumentError("target_faces must contain at least two entries"))

    for k in 1:length(source_faces)-1
        source_faces[k+1] > source_faces[k] ||
            throw(ArgumentError("source_faces must be strictly increasing"))
    end

    for k in 1:length(target_faces)-1
        target_faces[k+1] > target_faces[k] ||
            throw(ArgumentError("target_faces must be strictly increasing"))
    end

    return nothing
end

"""
    conservative_column_remap!(target_values, source_values, source_faces, target_faces; fill_value=0)

Conservatively remap piecewise-constant layer averages from `source_faces` to
`target_faces`. Uncovered parts of the target column are filled with
`fill_value`, which represents newly exposed or newly grown ice.
"""
function conservative_column_remap!(target_values,
                                    source_values,
                                    source_faces,
                                    target_faces;
                                    fill_value = zero(eltype(target_values)))
    check_column_remap_inputs(source_values, source_faces, target_faces)

    length(target_values) == length(target_faces) - 1 ||
        throw(ArgumentError("target_values must have one fewer entry than target_faces"))

    for kt in eachindex(target_values)
        target_bottom = target_faces[kt]
        target_top = target_faces[kt+1]
        target_width = target_top - target_bottom
        integral = fill_value * target_width

        for ks in eachindex(source_values)
            source_bottom = source_faces[ks]
            source_top = source_faces[ks+1]
            overlap = min(target_top, source_top) - max(target_bottom, source_bottom)

            if overlap > 0
                integral += (source_values[ks] - fill_value) * overlap
            end
        end

        target_values[kt] = integral / target_width
    end

    return target_values
end

"""
    conservative_column_remap(source_values, source_faces, target_faces; fill_value=0)

Return layer averages on `target_faces` by conservative piecewise-constant
remapping from `source_faces`.
"""
function conservative_column_remap(source_values,
                                   source_faces,
                                   target_faces;
                                   fill_value = zero(eltype(source_values)))
    target_values = similar(source_values, length(target_faces) - 1)
    return conservative_column_remap!(target_values,
                                      source_values,
                                      source_faces,
                                      target_faces;
                                      fill_value)
end

@inline column_remap_parameter(value::Number, i, j) = value
@inline column_remap_parameter(value::AbstractField, i, j) = @inbounds value[i, j, 1]

"""
    column_energy_thickness_remap!(thermodynamics, source_faces, target_faces; fill_energy=0, bulk_salinity=nothing)

Conservatively remap layer-averaged internal energy in `thermodynamics` from
`source_faces` to `target_faces`, then recompute thermodynamic diagnostics.
Uncovered target intervals are filled with `fill_energy`, which can be a scalar
or a two-dimensional field. If `bulk_salinity` is provided, the bulk-salinity
field is reset with that value or profile before diagnostics are recomputed;
this is the fixed-salinity path used after thickness changes.
"""
function column_energy_thickness_remap!(thermodynamics::ColumnEnergyThermodynamics,
                                        source_faces,
                                        target_faces;
                                        fill_energy = zero(eltype(thermodynamics.fields.internal_energy.grid)),
                                        bulk_salinity = nothing)
    grid = thermodynamics.fields.internal_energy.grid
    on_cpu(grid) || error("column_energy_thickness_remap! currently supports CPU grids only.")
    check_column_remap_inputs(zeros(eltype(grid), size(grid, 3)), source_faces, target_faces)

    Nz = size(grid, 3)
    length(target_faces) == Nz + 1 ||
        throw(ArgumentError("target_faces must have one more entry than the number of vertical cells"))

    fields = thermodynamics.fields
    source_values = zeros(eltype(grid), Nz)
    target_values = zeros(eltype(grid), Nz)

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        for k in 1:Nz
            @inbounds source_values[k] = fields.internal_energy[i, j, k]
        end

        fill_value = column_remap_parameter(fill_energy, i, j)
        conservative_column_remap!(target_values,
                                   source_values,
                                   source_faces,
                                   target_faces;
                                   fill_value)

        for k in 1:Nz
            @inbounds fields.internal_energy[i, j, k] = target_values[k]
        end
    end

    isnothing(bulk_salinity) || set_column_field!(fields.bulk_salinity, bulk_salinity)
    compute_column_thermodynamic_diagnostics!(thermodynamics)

    return nothing
end

"""
    column_layer_integral(values, faces)

Return the integral of piecewise-constant layer averages over `faces`.
"""
function column_layer_integral(values, faces)
    length(values) >= 1 ||
        throw(ArgumentError("values must contain at least one entry"))

    length(faces) == length(values) + 1 ||
        throw(ArgumentError("faces must have one more entry than values"))

    FT = promote_type(eltype(values), eltype(faces))
    integral = zero(FT)

    for k in eachindex(values)
        width = faces[k+1] - faces[k]
        width > 0 ||
            throw(ArgumentError("faces must be strictly increasing"))
        integral += values[k] * width
    end

    return integral
end

# Basal Stefan growth flux (positive = congelation). A `FluxBoundary` base injects its flux into the interior, so
# it drives no separate basal growth. A Dirichlet base (ocean/prescribed) grows ice from the conductive flux into
# the base minus the paired ocean heat flux Qᵇ — the same `(Q_cond − Q_ocean)/ℒ` balance SlabThermodynamics uses.
@inline column_basal_stefan_flux(boundary, ext, fields, auxiliary, grid, relation, clock, model_fields, i, j) = zero(eltype(grid))

@inline function column_basal_stefan_flux(boundary::ColumnDirichletBoundary, ext, fields, auxiliary, grid, relation, clock, model_fields, i, j)
    conductance = column_boundary_temperature_conductance(i, j, 1, 1, grid, auxiliary)
    T₁ = @inbounds fields.temperature[i, j, 1]
    Tᵇ = column_dirichlet_temperature(boundary, i, j, grid, relation)
    Qᵇ = getflux(ext, i, j, grid, Tᵇ, clock, model_fields)
    return conductance * (Tᵇ - T₁) - Qᵇ
end

# Unconsolidated top heat loss (upward-positive), mirroring SlabThermodynamics so a thin column behaves like a
# slab. A direct-flux/melting top exchanges its external flux; a Dirichlet (prescribed-temperature or ocean) top
# additionally loses heat by single-layer conduction k/h (Tᵇ − Tᵗ) — the same conduction the slab wraps into its
# top flux — so a thin conduction-driven column still grows instead of stalling.
@inline column_unconsolidated_top_flux(bc, ext, conductivity, Tᵇ, h, i, j, grid, fields, relation, clock, model_fields) =
    getflux(ext, i, j, grid, Tᵇ, clock, model_fields)

@inline function column_unconsolidated_top_flux(bc::ColumnDirichletBoundary, ext, conductivity, Tᵇ, h, i, j, grid, fields, relation, clock, model_fields)
    Tᵗ = column_dirichlet_temperature(bc, i, j, grid, relation)
    S  = @inbounds fields.bulk_salinity[i, j, 1]
    k  = ice_thermal_conductivity(conductivity, Tᵗ, S)
    Qi = ifelse(h ≤ zero(h), zero(h), k / h * (Tᵇ - Tᵗ))
    return Qi + getflux(ext, i, j, grid, Tᵗ, clock, model_fields)
end

@kernel function _column_stefan_volume_update!(ice_thickness, ice_concentration, consolidation_thickness,
                                               sea_ice_density, phase_transitions, z, fields, auxiliary, grid,
                                               heat_boundary_conditions, external_heat_fluxes, relation, conductivity,
                                               clock, model_fields, Δt)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        hⁿ = ice_thickness[i, j, 1]
        ℵⁿ = ice_concentration[i, j, 1]
        hᶜ = consolidation_thickness[i, j, 1]
        ρᵢ = column_stefan_parameter(sea_ice_density, i, j, 1)

        # Roll the interface heights into the previous slot for the next moving-grid solve.
        z.hs⁻[i, j, 1] = z.hsⁿ[i, j, 1]
        z.hb⁻[i, j, 1] = z.hbⁿ[i, j, 1]

        if hⁿ ≥ hᶜ
            # Consolidated: resolved conduction drives basal congelation and the surface melt residual.
            basal_flux = column_basal_stefan_flux(heat_boundary_conditions.bottom, external_heat_fluxes.bottom,
                                                  fields, auxiliary, grid, relation, clock, model_fields, i, j)
            surface_flux = column_surface_stefan_residual_flux(heat_boundary_conditions.top, external_heat_fluxes.top,
                                                               i, j, Nz, grid, auxiliary, fields, relation, clock, model_fields, Δt)

            δhᵇ = column_stefan_thickness_change(phase_transitions, ρᵢ, basal_flux, Δt)
            δhˢ = column_stefan_thickness_change(phase_transitions, ρᵢ, surface_flux, Δt)

            # Per-ice thickness velocity (no ℵ premultiplication), matching SlabThermodynamics — `ice_volume_update`
            # and the concentration step grow ℵ from open water, so a column nucleates ice from ℵ = 0 like the slab.
            ∂t_V = (δhᵇ + δhˢ) / Δt
            surface_shift = δhˢ
        else
            # Unconsolidated: no resolved interior profile, so the column behaves as a single-layer slab. Net growth
            # is ℰ⁻¹ (Qᵘ − Qᵇ) with the shared upward-positive fluxes, where a Dirichlet top contributes its
            # single-layer conduction k/h just as SlabThermodynamics does.
            Tᵇ = column_bottom_reference_temperature(heat_boundary_conditions.bottom, i, j, grid, fields, relation)
            ℰ = ρᵢ * latent_heat(phase_transitions, Tᵇ)
            Qᵘ = column_unconsolidated_top_flux(heat_boundary_conditions.top, external_heat_fluxes.top,
                                                conductivity, Tᵇ, hⁿ, i, j, grid, fields, relation, clock, model_fields)
            Qᵇ = getflux(external_heat_fluxes.bottom, i, j, grid, Tᵇ, clock, model_fields)

            ∂t_V = (Qᵘ - Qᵇ) / ℰ
            surface_shift = zero(eltype(grid))
        end

        h⁺, ℵ⁺ = ice_volume_update(ProportionalEvolution(), ∂t_V, hⁿ, ℵⁿ, hᶜ, Δt)

        # Keep the grid metric consistent with the (area-averaged) ice thickness: surface melt moves the surface,
        # and the base is placed so the column height equals the updated thickness. For ℵ = 1 this is exactly the
        # basal Stefan move (z.hbⁿ -= δhᵇ); for evolving ℵ it follows the redistributed thickness from ice_volume_update.
        # An ice-free column (h⁺ = 0) keeps its previous placeholder height so the metric never collapses to zero
        # (a zero Δz would make the resolved solve and the surface conductance singular).
        column_metric_height = ifelse(h⁺ > 0, h⁺, z.hs⁻[i, j, 1] - z.hb⁻[i, j, 1])
        z.hsⁿ[i, j, 1] = z.hs⁻[i, j, 1] + surface_shift
        z.hbⁿ[i, j, 1] = z.hsⁿ[i, j, 1] - column_metric_height

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
            thermal_conductivity(thermodynamics.energy_transport),
            model.clock, fields(model), Δt)
    return nothing
end

