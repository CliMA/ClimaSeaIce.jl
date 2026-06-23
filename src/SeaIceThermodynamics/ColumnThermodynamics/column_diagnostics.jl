#####
##### Diagnostics and column energy step
#####

thermal_conductivity(transport::ConductiveTemperatureTransport) = transport.conductivity
thermal_conductivity(::DiffusiveEnergyTransport) = 0
thermal_conductivity(transport::ConductiveAndDiffusiveEnergyTransport) = transport.conductivity

energy_diffusivity(::ConductiveTemperatureTransport) = 0
energy_diffusivity(transport::DiffusiveEnergyTransport) = transport.diffusivity
energy_diffusivity(transport::ConductiveAndDiffusiveEnergyTransport) = transport.diffusivity

salinity_diffusivity(::NoSalinityTransport) = 0
salinity_diffusivity(transport::BulkSalinityDiffusion) = transport.diffusivity

@inline shortwave_flux(::NoShortwaveAbsorption, z, z_top) = zero(z)

@inline shortwave_flux(absorption::ExponentialShortwaveAbsorption, z, z_top) =
    absorption.surface_transmission * exp((z - z_top) / absorption.attenuation_scale)

@kernel function _compute_column_internal_energy!(Φ, ℝ)
    i, j, k = @index(Global, NTuple)
    E = Φ.internal_energy
    
    @inbounds begin
        T = Φ.temperature[i, j, k]
        S = Φ.bulk_salinity[i, j, k]
        E[i, j, k] = internal_energy(ℝ, T, S)
    end
end

"""
    compute_column_internal_energy!(thermodynamics)

Fill the internal-energy field from the current temperature and bulk-salinity
fields.
"""
function compute_column_internal_energy!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_internal_energy!,
            thermodynamics.fields,
            thermodynamics.relation)

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

Fill temperature, liquid fraction, brine salinity, and thermodynamic derivative
diagnostics from internal energy and bulk salinity.
"""
function compute_column_thermodynamic_diagnostics!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_thermodynamic_diagnostics!,
            thermodynamics.fields,
            thermodynamics.relation)

    return nothing
end

@kernel function _compute_column_transport_coefficients!(auxiliary, fields, grid, energy_transport)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)

    κT = thermal_conductivity(energy_transport)
    κE = energy_diffusivity(energy_transport)

    @inbounds begin
        T = fields.temperature
        S = fields.bulk_salinity

        if k == 1
            boundary_conductivity = ice_thermal_conductivity(κT, T[i, j, 1], S[i, j, 1])
            TE = fields.temperature_energy_derivative[i, j, 1]
            TS = fields.temperature_salinity_derivative[i, j, 1]

            auxiliary.thermal_conductivity[i, j, 1] = boundary_conductivity
            auxiliary.energy_diffusivity[i, j, 1] = κE
            auxiliary.effective_energy_diffusivity[i, j, 1] = boundary_conductivity * TE + κE
            auxiliary.salinity_coupling_diffusivity[i, j, 1] = boundary_conductivity * TS
        end

        if k > 1
            TE = (fields.temperature_energy_derivative[i, j, k-1] +
                  fields.temperature_energy_derivative[i, j, k]) / 2

            TS = (fields.temperature_salinity_derivative[i, j, k-1] +
                  fields.temperature_salinity_derivative[i, j, k]) / 2

            face_conductivity =
                face_thermal_conductivity(κT,
                                          T[i, j, k-1],
                                          S[i, j, k-1],
                                          Δzᶜᶜᶜ(i, j, k-1, grid),
                                          T[i, j, k],
                                          S[i, j, k],
                                          Δzᶜᶜᶜ(i, j, k, grid))

            auxiliary.thermal_conductivity[i, j, k] = face_conductivity
            auxiliary.energy_diffusivity[i, j, k] = κE
            auxiliary.effective_energy_diffusivity[i, j, k] = face_conductivity * TE + κE
            auxiliary.salinity_coupling_diffusivity[i, j, k] = face_conductivity * TS
        end

        if k == Nz
            kp1 = Nz + 1
            boundary_conductivity = ice_thermal_conductivity(κT, T[i, j, Nz], S[i, j, Nz])
            TE = fields.temperature_energy_derivative[i, j, Nz]
            TS = fields.temperature_salinity_derivative[i, j, Nz]

            auxiliary.thermal_conductivity[i, j, kp1] = boundary_conductivity
            auxiliary.energy_diffusivity[i, j, kp1] = κE
            auxiliary.effective_energy_diffusivity[i, j, kp1] = boundary_conductivity * TE + κE
            auxiliary.salinity_coupling_diffusivity[i, j, kp1] = boundary_conductivity * TS
        end
    end
end

"""
    compute_column_transport_coefficients!(thermodynamics)

Compute face-centered energy transport coefficients for the semi-implicit
energy solve.
"""
function compute_column_transport_coefficients!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_transport_coefficients!,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            grid,
            thermodynamics.energy_transport)

    return nothing
end

@kernel function _compute_column_salinity_diffusivity!(auxiliary, fields, salinity_transport)
    i, j, k = @index(Global, NTuple)
    Nz = size(fields.bulk_salinity.grid, 3)

    κS = salinity_diffusivity(salinity_transport)

    @inbounds begin
        if k == 1
            auxiliary.salinity_diffusivity[i, j, 1] = zero(κS)
        end

        if k > 1
            auxiliary.salinity_diffusivity[i, j, k] = κS
        end

        if k == Nz
            auxiliary.salinity_diffusivity[i, j, Nz+1] = zero(κS)
        end
    end
end

"""
    compute_column_salinity_diffusivity!(thermodynamics)

Compute face-centered bulk-salinity diffusivity for the scalar salinity solve.
"""
function compute_column_salinity_diffusivity!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.bulk_salinity.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_salinity_diffusivity!,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            thermodynamics.salinity_transport)

    return nothing
end

@kernel function _compute_column_shortwave_flux!(auxiliary, grid, shortwave_absorption)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)
    z_top = znode(i, j, Nz+1, grid, Center(), Center(), Face())

    @inbounds begin
        if k == 1
            z_bottom = znode(i, j, 1, grid, Center(), Center(), Face())
            auxiliary.shortwave_flux[i, j, 1] =
                shortwave_flux(shortwave_absorption, z_bottom, z_top)
        end

        if k > 1
            z_face = znode(i, j, k, grid, Center(), Center(), Face())
            auxiliary.shortwave_flux[i, j, k] =
                shortwave_flux(shortwave_absorption, z_face, z_top)
        end

        if k == Nz
            auxiliary.shortwave_flux[i, j, Nz+1] =
                shortwave_flux(shortwave_absorption, z_top, z_top)
        end
    end
end

"""
    compute_column_shortwave_flux!(thermodynamics)

Compute face-centered shortwave fluxes for the current shortwave absorption
closure.
"""
function compute_column_shortwave_flux!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_shortwave_flux!,
            thermodynamics.auxiliary,
            grid,
            thermodynamics.shortwave_absorption)

    return nothing
end

@kernel function _compute_column_surface_stefan_residual_flux!(residual_flux,
                                                               auxiliary,
                                                               fields,
                                                               grid,
                                                               heat_boundary_conditions,
                                                               external_heat_fluxes,
                                                               relation,
                                                               clock,
                                                               model_fields,
                                                               Δt)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds residual_flux[i, j, 1] =
        column_surface_stefan_residual_flux(heat_boundary_conditions.top, external_heat_fluxes.top,
                                            i, j, Nz, grid, auxiliary, fields, relation, clock, model_fields, Δt)
end

"""
    compute_column_surface_stefan_residual_flux!(residual_flux, thermodynamics, external_heat_fluxes, clock, model_fields, dt)

Fill `residual_flux` with the top-surface Stefan residual implied by a
`MeltingConstrainedFluxBalance` surface boundary. The residual is positive for
growth and negative for melt, matching [`column_stefan_thickness_update!`](@ref).
"""
function compute_column_surface_stefan_residual_flux!(residual_flux,
                                                      thermodynamics::ColumnEnergyThermodynamics,
                                                      external_heat_fluxes,
                                                      clock,
                                                      model_fields,
                                                      Δt)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    # The residual caps against the start-of-step (current) top-cell enthalpy; stash it so the residual reads the
    # same enthalpy whether called standalone here or after the coupled surface solve has advanced internal_energy.
    Nz = size(grid, 3)
    interior(thermodynamics.auxiliary.surface_start_energy)[:, :, 1] .=
        @view interior(thermodynamics.fields.internal_energy)[:, :, Nz]

    launch!(arch, grid, :xy,
            _compute_column_surface_stefan_residual_flux!,
            residual_flux,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            grid,
            thermodynamics.heat_boundary_conditions,
            external_heat_fluxes,
            thermodynamics.relation,
            clock,
            model_fields,
            Δt)

    return nothing
end

@inline diffusion_factor(i, j, kᶜ, kᶠ, grid, K, Δt) =  Δt * @inbounds(K[i, j, kᶠ]) / (Δzᶜᶜᶜ(i, j, kᶜ, grid) * Δzᶜᶜᶠ(i, j, kᶠ, grid))

# Conservative metric Jacobian: the column-thickness ratio (σ ≡ column height / reference height).
@inline previous_to_current_column_metric_ratio(i, j, k, grid::SeaIceColumnGrid) =
    previous_column_height(grid, i, j) / column_height(grid, i, j)

# The two interfaces move independently, so the swept-face displacement is the base motion plus the standard
# thickness-stretch term (σⁿ - σ⁻) r: basal growth/melt moves only the base, surface melt/snow-ice only the surface.
@inline function column_face_displacement(i, j, k, grid::SeaIceColumnGrid)
    z = sea_ice_discretization(grid)
    δbase = @inbounds z.hbⁿ[i, j, 1] - z.hb⁻[i, j, 1]
    σᶜᶜⁿ = σⁿ(i, j, k, grid, Center(), Center(), Face())
    σᶜᶜ⁻ = σ⁻(i, j, k, grid, Center(), Center(), Face())
    r = rnode(i, j, k, grid, Center(), Center(), Face())
    return δbase + (σᶜᶜⁿ - σᶜᶜ⁻) * r
end

# Dirichlet boundaries (`PrescribedTemperature`, `IceWaterThermalEquilibrium`) add an implicit conductance factor
# and a conductive flux; every other boundary injects its paired external flux directly.
const ColumnDirichletBoundary = Union{PrescribedTemperature, IceWaterThermalEquilibrium}

# The temperature of a Dirichlet boundary (ocean equilibrium or prescribed), shared with the slab vocabulary.
@inline column_dirichlet_temperature(bc, i, j, grid, relation) =
    bottom_temperature(i, j, grid, bc, relation.phase_transitions.liquidus)

# A representative bottom temperature for the latent-heat scaling: the Dirichlet interface temperature where one
# exists, otherwise the adjacent bottom-cell temperature (e.g. a `FluxBoundary` base driven by a prescribed flux).
@inline column_bottom_reference_temperature(bc, i, j, grid, fields, relation) = @inbounds fields.temperature[i, j, 1]
@inline column_bottom_reference_temperature(bc::ColumnDirichletBoundary, i, j, grid, fields, relation) =
    column_dirichlet_temperature(bc, i, j, grid, relation)

# Swept material on a moving interface takes the adjacent interior enthalpy by default. A Dirichlet basal interface
# (congelation) instead lays down ice at the interface temperature — the CICE BL99 new-ice enthalpy.
@inline moving_bottom_boundary_value(boundary, i, j, grid, fields, relation, default) = default
@inline moving_top_boundary_value(boundary, i, j, grid, fields, relation, default) = default

@inline function moving_bottom_boundary_value(boundary::Union{PrescribedTemperature, IceWaterThermalEquilibrium}, i, j, grid, fields, relation, default)
    Tᵇ = column_dirichlet_temperature(boundary, i, j, grid, relation)
    S = @inbounds fields.bulk_salinity[i, j, 1]
    return internal_energy(relation, Tᵇ, S)
end

@inline function moving_face_value(i, j, k, grid, field, face_displacement, bottom_boundary_value, top_boundary_value)
    Nz = size(grid, 3)

    if k == 1
        interior_value = @inbounds field[i, j, 1]
        return face_displacement < 0 ? bottom_boundary_value : interior_value
    elseif k == Nz + 1
        interior_value = @inbounds field[i, j, Nz]
        return face_displacement > 0 ? top_boundary_value : interior_value
    elseif face_displacement >= 0
        return @inbounds field[i, j, k]
    else
        return @inbounds field[i, j, k-1]
    end
end

@inline function moving_face_displacement_flux(i, j, k, grid, field)
    displacement = column_face_displacement(i, j, k, grid)
    Nz = size(grid, 3)
    bottom_boundary_value = @inbounds field[i, j, 1]
    top_boundary_value = @inbounds field[i, j, Nz]
    value = moving_face_value(i, j, k, grid, field, displacement,
                              bottom_boundary_value, top_boundary_value)
    return displacement * value
end

@inline function moving_face_displacement_flux(i, j, k, grid, field, fields, relation,
                                               bottom_boundary,
                                               top_boundary)
    displacement = column_face_displacement(i, j, k, grid)
    Nz = size(grid, 3)
    bottom_default = @inbounds field[i, j, 1]
    top_default = @inbounds field[i, j, Nz]
    bottom_boundary_value = moving_bottom_boundary_value(bottom_boundary, i, j, grid, fields, relation, bottom_default)
    top_boundary_value = moving_top_boundary_value(top_boundary, i, j, grid, fields, relation, top_default)
    value = moving_face_value(i, j, k, grid, field, displacement, bottom_boundary_value, top_boundary_value)
    return displacement * value
end

@inline salinity_coupling_flux(i, j, k, grid, C, S) = @inbounds(C[i, j, k] * (S[i, j, k] - S[i, j, k-1]) / Δzᶜᶜᶠ(i, j, k, grid))

@inline column_boundary_temperature_conductance(i, j, kf, kc, grid, auxiliary) = @inbounds 2 * auxiliary.thermal_conductivity[i, j, kf] / Δzᶜᶜᶜ(i, j, kc, grid)

@inline column_consolidated(consolidation_thickness::Nothing, grid, i, j) = true
@inline function column_consolidated(consolidation_thickness, grid, i, j)
    h = column_height(grid, i, j)
    hᶜ = @inbounds consolidation_thickness[i, j, 1]
    return h ≥ hᶜ
end

@inline function column_energy_to_complete_melt(relation, E, S, Δz)
    return (complete_melt_energy(relation, S) - E) * Δz
end

