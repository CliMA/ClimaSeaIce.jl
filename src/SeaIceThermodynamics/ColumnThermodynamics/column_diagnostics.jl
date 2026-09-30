#####
##### Diagnostics and column energy step
#####

thermal_conductivity(transport::ConductiveTemperatureTransport) = transport.conductivity
energy_diffusivity(transport::ConductiveTemperatureTransport) = transport.diffusivity

salinity_diffusivity(::NoSalinityTransport) = 0
salinity_diffusivity(transport::BulkSalinityDiffusion) = transport.diffusivity

@inline shortwave_flux(::NoShortwaveAbsorption, z, z_top) = zero(z)
@inline shortwave_flux(absorption::ExponentialShortwaveAbsorption, z, z_top) = absorption.surface_transmission * exp((z - z_top) / absorption.attenuation_scale)

@kernel function _compute_column_internal_energy!(Φ, ℝ)
    i, j, k = @index(Global, NTuple)
    @inbounds Φ.internal_energy[i, j, k] = internal_energy(ℝ, Φ.temperature[i, j, k], Φ.bulk_salinity[i, j, k])
end

"""
    compute_column_internal_energy!(thermodynamics)

Fill the internal-energy field from the current temperature and bulk-salinity fields.
"""
function compute_column_internal_energy!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xyz, _compute_column_internal_energy!, thermodynamics.fields, thermodynamics.relation)
    return nothing
end

@kernel function _compute_column_thermodynamic_diagnostics!(Φ, ℝ)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        E = Φ.internal_energy[i, j, k]
        S = Φ.bulk_salinity[i, j, k]

        Φ.temperature[i, j, k] = temperature(ℝ, E, S)
        Φ.liquid_fraction[i, j, k] = liquid_fraction(ℝ, E, S)
        Φ.brine_salinity[i, j, k] = brine_salinity(ℝ, E, S)
        Φ.temperature_energy_derivative[i, j, k] = temperature_energy_derivative(ℝ, E, S)
        Φ.temperature_salinity_derivative[i, j, k] = temperature_salinity_derivative(ℝ, E, S)
    end
end

"""
    compute_column_thermodynamic_diagnostics!(thermodynamics)

Fill temperature, liquid fraction, brine salinity, and the thermodynamic derivatives from internal energy and bulk salinity.
"""
function compute_column_thermodynamic_diagnostics!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xyz, _compute_column_thermodynamic_diagnostics!, thermodynamics.fields, thermodynamics.relation)
    return nothing
end

# Boundary faces (k = 1 and k = Nz + 1) take the adjacent cell values.
@inline function store_face_transport_coefficients!(auxiliary, i, j, k, grid, fields, κT, κE)
    Nz = size(grid, 3)
    kᵇ = max(k - 1, 1)
    kᵗ = min(k, Nz)

    @inbounds begin
        T = fields.temperature
        S = fields.bulk_salinity
        TE = (fields.temperature_energy_derivative[i, j, kᵇ] + fields.temperature_energy_derivative[i, j, kᵗ]) / 2
        TS = (fields.temperature_salinity_derivative[i, j, kᵇ] + fields.temperature_salinity_derivative[i, j, kᵗ]) / 2

        K = ifelse(kᵇ == kᵗ, ice_thermal_conductivity(κT, T[i, j, kᵗ], S[i, j, kᵗ]),
                   face_thermal_conductivity(κT, T[i, j, kᵇ], S[i, j, kᵇ], Δzᶜᶜᶜ(i, j, kᵇ, grid), T[i, j, kᵗ], S[i, j, kᵗ], Δzᶜᶜᶜ(i, j, kᵗ, grid)))

        auxiliary.thermal_conductivity[i, j, k] = K
        auxiliary.effective_energy_diffusivity[i, j, k] = K * TE + κE
        auxiliary.salinity_coupling_diffusivity[i, j, k] = K * TS
    end
end

@kernel function _compute_column_transport_coefficients!(auxiliary, fields, grid, energy_transport)
    i, j, k = @index(Global, NTuple)
    κT = thermal_conductivity(energy_transport)
    κE = energy_diffusivity(energy_transport)
    store_face_transport_coefficients!(auxiliary, i, j, k, grid, fields, κT, κE)
    k == size(grid, 3) && store_face_transport_coefficients!(auxiliary, i, j, k + 1, grid, fields, κT, κE)
end

"""
    compute_column_transport_coefficients!(thermodynamics)

Compute the face-centered energy transport coefficients of the semi-implicit energy solve.
"""
function compute_column_transport_coefficients!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xyz, _compute_column_transport_coefficients!,
            thermodynamics.auxiliary, thermodynamics.fields, grid, thermodynamics.energy_transport)
    return nothing
end

@kernel function _compute_column_shortwave_flux!(I, grid, shortwave_absorption)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)
    z_top = znode(i, j, Nz+1, grid, Center(), Center(), Face())
    @inbounds I[i, j, k] = shortwave_flux(shortwave_absorption, znode(i, j, k, grid, Center(), Center(), Face()), z_top)
    @inbounds k == Nz && (I[i, j, Nz+1] = shortwave_flux(shortwave_absorption, z_top, z_top))
end

"""
    compute_column_shortwave_flux!(thermodynamics)

Compute the face-centered shortwave flux of the shortwave absorption closure.
"""
function compute_column_shortwave_flux!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xyz, _compute_column_shortwave_flux!,
            thermodynamics.auxiliary.shortwave_flux, grid, thermodynamics.shortwave_absorption)
    return nothing
end

@inline diffusion_factor(i, j, kᶜ, kᶠ, grid, K::Number, Δt) = Δt * K / (Δzᶜᶜᶜ(i, j, kᶜ, grid) * Δzᶜᶜᶠ(i, j, kᶠ, grid))
@inline diffusion_factor(i, j, kᶜ, kᶠ, grid, K, Δt) = diffusion_factor(i, j, kᶜ, kᶠ, grid, @inbounds(K[i, j, kᶠ]), Δt)

@inline previous_to_current_column_metric_ratio(i, j, k, grid::SeaIceColumnGrid) = previous_column_height(grid, i, j) / column_height(grid, i, j)

# Swept-face displacement: base motion plus the thickness-stretch term (σⁿ - σ⁻) r.
@inline function column_face_displacement(i, j, k, grid::SeaIceColumnGrid)
    z = sea_ice_discretization(grid)
    δbase = @inbounds z.hbⁿ[i, j, 1] - z.hb⁻[i, j, 1]
    σᶜᶜⁿ = σⁿ(i, j, k, grid, Center(), Center(), Face())
    σᶜᶜ⁻ = σ⁻(i, j, k, grid, Center(), Center(), Face())
    return δbase + (σᶜᶜⁿ - σᶜᶜ⁻) * rnode(i, j, k, grid, Center(), Center(), Face())
end

const ColumnDirichletBoundary = Union{PrescribedTemperature, IceWaterThermalEquilibrium}

@inline column_dirichlet_temperature(bc, i, j, grid, relation) = bottom_temperature(i, j, grid, bc, relation.phase_transitions.liquidus)

@inline column_bottom_reference_temperature(bc, i, j, grid, fields, relation) = @inbounds fields.temperature[i, j, 1]
@inline column_bottom_reference_temperature(bc::ColumnDirichletBoundary, i, j, grid, fields, relation) = column_dirichlet_temperature(bc, i, j, grid, relation)

# Congelation at a Dirichlet base lays down ice at the interface temperature (the BL99 new-ice enthalpy).
@inline moving_bottom_boundary_value(boundary, i, j, grid, fields, relation, default) = default

@inline function moving_bottom_boundary_value(boundary::ColumnDirichletBoundary, i, j, grid, fields, relation, default)
    Tᵇ = column_dirichlet_temperature(boundary, i, j, grid, relation)
    return internal_energy(relation, Tᵇ, @inbounds(fields.bulk_salinity[i, j, 1]))
end

# Upwind value swept across face k; boundary faces sweep in the bottom (top) boundary value when the ice grows.
@inline function moving_face_displacement_flux(i, j, k, grid, c, fields, relation, bottom_boundary)
    Nz = size(grid, 3)
    δz = column_face_displacement(i, j, k, grid)

    @inbounds begin
        cᵇ = moving_bottom_boundary_value(bottom_boundary, i, j, grid, fields, relation, c[i, j, 1])
        interior_value = ifelse(δz >= 0, c[i, j, min(k, Nz)], c[i, j, max(k-1, 1)])
        value = ifelse((k == 1) & (δz < 0), cᵇ, interior_value)
    end

    return δz * value
end

@inline salinity_coupling_flux(i, j, k, grid, C, S) = @inbounds(C[i, j, k] * (S[i, j, k] - S[i, j, k-1]) / Δzᶜᶜᶠ(i, j, k, grid))

@inline column_boundary_temperature_conductance(i, j, kf, kc, grid, auxiliary) = @inbounds 2 * auxiliary.thermal_conductivity[i, j, kf] / Δzᶜᶜᶜ(i, j, kc, grid)

@inline column_consolidated(consolidation_thickness, grid, i, j) = column_height(grid, i, j) ≥ @inbounds consolidation_thickness[i, j, 1]

@inline column_energy_to_complete_melt(relation, E, S, Δz) = (complete_melt_energy(relation, S) - E) * Δz
