using ClimaSeaIce.SeaIceThermodynamics:
    BulkSalinityDiffusion,
    ColumnEnergyThermodynamics,
    SeaIceColumnDiscretization,
    ConductiveTemperatureTransport,
    ExponentialShortwaveAbsorption,
    BubblyBrineConductivity,
    FixedDrainedIceSalinityProfile,
    FixedSalinityBrinePocketEnergyRelation,
    FluxBoundary,
    FluxFunction,
    MeltingConstrainedFluxBalance,
    IceWaterThermalEquilibrium,
    MaykutUntersteinerConductivity,
    NoSalinityTransport,
    NoShortwaveAbsorption,
    PrescribedBulkSalinity,
    PrescribedTemperature,
    PrognosticBulkSalinity,
    QuadraticLiquidusEnergyRelation,
    assemble_column_energy_system!,
    brine_salinity,
    column_energy_time_step!,
    column_height,
    column_salinity_time_step!,
    column_surface_stefan_residual_flux,
    column_stefan_thickness_change,
    complete_melt_energy,
    compute_column_shortwave_flux!,
    compute_column_thermodynamic_diagnostics!,
    compute_column_transport_coefficients!,
    evolving_salinity_mushy_thermodynamics,
    face_thermal_conductivity,
    ice_thermal_conductivity,
    internal_energy,
    liquid_fraction,
    prescribed_salinity_enthalpy_thermodynamics,
    salinity_at_normalized_depth,
    salinity_at_normalized_height,
    solve_column_energy_system!,
    temperature,
    temperature_energy_derivative,
    temperature_salinity_derivative

using ClimaSeaIce: SeaIceModel
using Adapt
using Oceananigans
using Oceananigans.TimeSteppers: Clock, time_step!
using Oceananigans.Fields: interior, set!, ConstantField
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans: fields, prognostic_fields, prognostic_state, restore_prognostic_state!
using Test

# Unit tests drive the column step with a model-style `external_heat_fluxes = (top, bottom)` set, a clock, no model
# fields, and a zero consolidation thickness so that every column is consolidated.
const TEST_CLOCK = Clock(time = 0.0)
const NO_FIELDS = NamedTuple()
const CONSOLIDATED = ConstantField(0.0)

insulating_fluxes() = (top = 0.0, bottom = 0.0)
insulating_boundaries() = (top = FluxBoundary(), bottom = FluxBoundary())

column_step!(thermodynamics, external_heat_fluxes, Δt) =
    column_energy_time_step!(thermodynamics, external_heat_fluxes, TEST_CLOCK, NO_FIELDS, CONSOLIDATED, Δt)

column_integral(c) = sum(interior(Field(c * KernelFunctionOperation{Center, Center, Center}(Δzᶜᶜᶜ, c.grid))))
column_integrated_energy(thermodynamics) = column_integral(thermodynamics.fields.internal_energy)
column_integrated_salinity(thermodynamics) = column_integral(thermodynamics.fields.bulk_salinity)

# Reference conservative remap of piecewise-constant layer averages; uncovered target intervals take `fill_value`.
function conservative_column_remap(source_values, source_faces, target_faces; fill_value = zero(eltype(source_values)))
    return map(1:length(target_faces)-1) do k
        zᵇ, zᵗ = target_faces[k], target_faces[k+1]
        overlaps = [max(min(zᵗ, source_faces[s+1]) - max(zᵇ, source_faces[s]), 0) for s in eachindex(source_values)]
        (fill_value * (zᵗ - zᵇ) + sum((source_values .- fill_value) .* overlaps)) / (zᵗ - zᵇ)
    end
end

function column_relation_test_states(::Type{FT}) where FT
    temperatures = FT[-40, -20, -5, -1, -0.1]
    salinities = FT[0, 0.5, 4, 8, 15]
    return temperatures, salinities
end

function physical_column_relation_test_states(::Type{FT}) where FT
    salinities = FT[0.1, 1, 4, 8, 15]
    offsets = FT[0.05, 0.5, 2, 10, 30]
    return salinities, offsets
end

function test_temperature_inversion(::Type{FT}) where FT
    relation = QuadraticLiquidusEnergyRelation(FT)
    temperatures, salinities = column_relation_test_states(FT)

    max_error = zero(FT)
    for T in temperatures, S in salinities
        E = internal_energy(relation, T, S)
        recovered_T = temperature(relation, E, S)
        max_error = max(max_error, abs(recovered_T - T))
    end

    return max_error
end

function test_temperature_derivatives(::Type{FT}) where FT
    relation = QuadraticLiquidusEnergyRelation(FT)
    WT = FT === Float32 ? Float64 : FT
    T0 = relation.phase_transitions.reference_temperature
    m = relation.phase_transitions.liquidus.salinity_slope
    salinities, offsets = physical_column_relation_test_states(FT)

    max_energy_derivative_error = zero(FT)
    max_salinity_derivative_error = zero(FT)

    for S in salinities, offset in offsets
        T = T0 - m * S - offset
        E = internal_energy(relation, T, S)

        dE = sqrt(eps(WT)) * max(abs(WT(E)), one(WT))
        dS = sqrt(eps(WT)) * max(abs(WT(S)), one(WT))

        analytic_TE = temperature_energy_derivative(relation, E, S)
        fd_TE = (temperature(relation, WT(E) + dE, WT(S)) -
                 temperature(relation, WT(E) - dE, WT(S))) / (2 * dE)

        analytic_TS = temperature_salinity_derivative(relation, E, S)
        fd_TS = (temperature(relation, WT(E), WT(S) + dS) -
                 temperature(relation, WT(E), WT(S) - dS)) / (2 * dS)

        TE_error = abs(analytic_TE - fd_TE) /
                   max(abs(analytic_TE), abs(fd_TE), one(FT))

        TS_error = abs(analytic_TS - fd_TS) /
                   max(abs(analytic_TS), abs(fd_TS), one(FT))

        max_energy_derivative_error = max(max_energy_derivative_error, TE_error)
        max_salinity_derivative_error = max(max_salinity_derivative_error, TS_error)
    end

    return max(max_energy_derivative_error, max_salinity_derivative_error)
end

function test_pure_ice_limit(::Type{FT}) where FT
    relation = QuadraticLiquidusEnergyRelation(FT)
    phase_transitions = relation.phase_transitions
    T0 = phase_transitions.reference_temperature
    heat_capacity = phase_transitions.density * phase_transitions.heat_capacity
    temperatures = FT[-40, -20, -5, -1, -0.1, 0]
    S = zero(FT)

    max_temperature_error = zero(FT)
    max_energy_error = zero(FT)
    max_derivative_error = zero(FT)

    for T in temperatures
        expected_E = heat_capacity * (T - T0)
        E = internal_energy(relation, T, S)
        recovered_T = temperature(relation, E, S)
        TE = temperature_energy_derivative(relation, E, S)

        max_temperature_error = max(max_temperature_error, abs(recovered_T - T))
        max_energy_error = max(max_energy_error, abs(E - expected_E))
        max_derivative_error = max(max_derivative_error, abs(TE - inv(heat_capacity)))
    end

    return max_temperature_error, max_energy_error, max_derivative_error
end

function test_liquidus_consistency(::Type{FT}) where FT
    relation = QuadraticLiquidusEnergyRelation(FT)
    liquidus = relation.phase_transitions.liquidus
    T0 = relation.phase_transitions.reference_temperature
    m = liquidus.salinity_slope
    salinities, offsets = physical_column_relation_test_states(FT)

    max_error = zero(FT)
    for S in salinities, offset in offsets
        T = T0 - m * S - offset
        E = internal_energy(relation, T, S)
        brine_S = brine_salinity(relation, E, S)
        recovered_T = temperature(relation, E, S)
        liquidus_T = Float64(liquidus.freshwater_melting_temperature) -
                     Float64(liquidus.salinity_slope) * Float64(brine_S)
        max_error = max(max_error, abs(Float64(recovered_T) - liquidus_T))
    end

    return max_error
end

function test_liquid_fraction_bounds(::Type{FT}) where FT
    relation = QuadraticLiquidusEnergyRelation(FT)
    T0 = relation.phase_transitions.reference_temperature
    m = relation.phase_transitions.liquidus.salinity_slope
    salinities, offsets = physical_column_relation_test_states(FT)

    for S in salinities, offset in offsets
        T = T0 - m * S - offset
        E = internal_energy(relation, T, S)
        phi = liquid_fraction(relation, E, S)
        brine_S = brine_salinity(relation, E, S)

        @test isfinite(phi)
        @test isfinite(brine_S)
        @test zero(FT) <= phi <= one(FT)
    end

    return nothing
end

function test_complete_melt_energy(::Type{FT}) where FT
    relation = QuadraticLiquidusEnergyRelation(FT)
    T0 = relation.phase_transitions.reference_temperature
    m = relation.phase_transitions.liquidus.salinity_slope
    salinities = FT[0.1, 1, 4, 8, 15]

    max_deficit = zero(FT)
    for S in salinities
        T = T0 - m * S
        E = internal_energy(relation, T, S)
        melt_E = complete_melt_energy(relation, S)
        scale = relation.phase_transitions.liquid_density *
                relation.phase_transitions.reference_latent_heat
        deficit = (melt_E - E) / scale
        max_deficit = max(max_deficit, abs(deficit))
    end

    return max_deficit
end

column_values(field) = vec(Array(interior(field)))

# Set the SeaIceColumnDiscretization interfaces to a previous→current thickness change. `move = :top` keeps the
# base at `surface` and moves the surface (the default, e.g. surface melt/snow-ice); `move = :base` keeps the
# surface at `surface` and moves the base (basal congelation/melt).
function set_column_metric!(grid, previous, current; surface = 0, move = :top)
    z = grid.z
    if move === :top
        fill!(z.hb⁻, surface); fill!(z.hbⁿ, surface)
        fill!(z.hs⁻, surface + previous); fill!(z.hsⁿ, surface + current)
    else
        fill!(z.hs⁻, surface); fill!(z.hsⁿ, surface)
        fill!(z.hb⁻, surface - previous); fill!(z.hbⁿ, surface - current)
    end
    return nothing
end

function bl99_salinity_reference(depth)
    return setprecision(BigFloat, 256) do
        s = BigFloat(depth)
        Smax = BigFloat("3.2")
        a = BigFloat("0.407")
        b = BigFloat("0.573")
        Smax / 2 * (1 - cos(big(pi) * s^(a / (s + b))))
    end
end

function test_fixed_salinity_temperature_inversion(::Type{FT}) where FT
    relation = FixedSalinityBrinePocketEnergyRelation(FT)
    temperatures, salinities = column_relation_test_states(FT)

    max_error = zero(FT)
    for T in temperatures, S in salinities
        E = internal_energy(relation, T, S)
        recovered_T = temperature(relation, E, S)
        max_error = max(max_error, abs(recovered_T - T))
    end

    return max_error
end

function test_fixed_salinity_temperature_derivatives(::Type{FT}) where FT
    relation = FixedSalinityBrinePocketEnergyRelation(FT)
    WT = FT === Float32 ? Float64 : FT
    T0 = relation.phase_transitions.reference_temperature
    m = relation.phase_transitions.liquidus.salinity_slope
    salinities, offsets = physical_column_relation_test_states(FT)

    max_energy_derivative_error = zero(FT)
    max_salinity_derivative_error = zero(FT)

    for S in salinities, offset in offsets
        T = T0 - m * S - offset
        E = internal_energy(relation, T, S)

        dE = sqrt(eps(WT)) * max(abs(WT(E)), one(WT))
        dS = sqrt(eps(WT)) * max(abs(WT(S)), one(WT))

        analytic_TE = temperature_energy_derivative(relation, E, S)
        fd_TE = (temperature(relation, WT(E) + dE, WT(S)) -
                 temperature(relation, WT(E) - dE, WT(S))) / (2 * dE)

        analytic_TS = temperature_salinity_derivative(relation, E, S)
        fd_TS = (temperature(relation, WT(E), WT(S) + dS) -
                 temperature(relation, WT(E), WT(S) - dS)) / (2 * dS)

        TE_error = abs(analytic_TE - fd_TE) /
                   max(abs(analytic_TE), abs(fd_TE), one(FT))

        TS_error = abs(analytic_TS - fd_TS) /
                   max(abs(analytic_TS), abs(fd_TS), one(FT))

        max_energy_derivative_error = max(max_energy_derivative_error, TE_error)
        max_salinity_derivative_error = max(max_salinity_derivative_error, TS_error)
    end

    return max(max_energy_derivative_error, max_salinity_derivative_error)
end

function test_fixed_salinity_pure_ice_limit(::Type{FT}) where FT
    relation = FixedSalinityBrinePocketEnergyRelation(FT)
    phase_transitions = relation.phase_transitions
    T0 = phase_transitions.reference_temperature
    heat_capacity = phase_transitions.density * phase_transitions.heat_capacity
    temperatures = FT[-40, -20, -5, -1, -0.1, 0]
    S = zero(FT)

    max_temperature_error = zero(FT)
    max_energy_error = zero(FT)
    max_derivative_error = zero(FT)

    for T in temperatures
        expected_E = heat_capacity * (T - T0)
        E = internal_energy(relation, T, S)
        recovered_T = temperature(relation, E, S)
        TE = temperature_energy_derivative(relation, E, S)

        max_temperature_error = max(max_temperature_error, abs(recovered_T - T))
        max_energy_error = max(max_energy_error, abs(E - expected_E))
        max_derivative_error = max(max_derivative_error, abs(TE - inv(heat_capacity)))
    end

    return max_temperature_error, max_energy_error, max_derivative_error
end

function test_fixed_salinity_complete_melt_threshold(::Type{FT}) where FT
    relation = FixedSalinityBrinePocketEnergyRelation(FT)
    T0 = relation.phase_transitions.reference_temperature
    m = relation.phase_transitions.liquidus.salinity_slope
    salinities = FT[0.1, 1, 3.2, 8, 15]

    max_error = zero(FT)
    for S in salinities
        melt_E = complete_melt_energy(relation, S)
        melt_T = temperature(relation, melt_E, S)
        max_error = max(max_error, abs(melt_T - (T0 - m * S)))
    end

    return max_error
end

function test_fixed_salinity_heat_capacity_identity(::Type{FT}) where FT
    relation = FixedSalinityBrinePocketEnergyRelation(FT)
    phase_transitions = relation.phase_transitions
    T0 = phase_transitions.reference_temperature
    m = phase_transitions.liquidus.salinity_slope
    rho_i = phase_transitions.density
    c_i = phase_transitions.heat_capacity
    L0 = phase_transitions.reference_latent_heat
    salinities, offsets = physical_column_relation_test_states(FT)

    max_error = zero(FT)
    for S in salinities, offset in offsets
        T = T0 - m * S - offset
        E = internal_energy(relation, T, S)
        tau = T - T0
        denominator = inv(temperature_energy_derivative(relation, E, S))
        expected = rho_i * (c_i + L0 * m * S / tau^2)
        error = abs(denominator - expected) / max(abs(denominator), abs(expected), one(FT))
        max_error = max(max_error, error)
    end

    return max_error
end

function maykut_untersteiner_reference(conductivity, temperature, salinity)
    T = min(-conductivity.temperature_floor, temperature)
    K = conductivity.fresh_ice_conductivity +
        conductivity.salinity_coefficient * salinity / T
    return max(K, conductivity.minimum_conductivity)
end

function bubbly_brine_reference(conductivity, temperature, salinity)
    T = min(-conductivity.temperature_floor, temperature)
    K = conductivity.ice_density / conductivity.pure_ice_density *
        (2.11 - 0.011 * temperature + 0.09 * salinity / T)
    return max(K, conductivity.minimum_conductivity)
end

function solve_reference_vector_tridiagonal(lower_diagonal, diagonal, upper_diagonal, rhs)
    n = length(rhs)
    solution = similar(rhs)
    scratch = zeros(eltype(rhs), max(n - 1, 0))
    β = first(diagonal)

    solution[1] = rhs[1] / β

    for k in 2:n
        scratch[k-1] = upper_diagonal[k-1] / β
        β = diagonal[k] - lower_diagonal[k-1] * scratch[k-1]
        solution[k] = (rhs[k] - lower_diagonal[k-1] * solution[k-1]) / β
    end

    for k in n-1:-1:1
        solution[k] -= scratch[k] * solution[k+1]
    end

    return solution
end

function dense_column_energy_matrix(thermodynamics)
    lower = column_values(thermodynamics.auxiliary.lower_diagonal)
    diagonal = column_values(thermodynamics.auxiliary.diagonal)
    upper = column_values(thermodynamics.auxiliary.upper_diagonal)
    Nz = length(diagonal)
    matrix = zeros(eltype(diagonal), Nz, Nz)

    for k in 1:Nz
        matrix[k, k] = diagonal[k]

        if k > 1
            matrix[k, k-1] = lower[k-1]
        end

        if k < Nz
            matrix[k, k+1] = upper[k]
        end
    end

    return matrix
end

function pure_ice_conductive_column(N; diffusivity = 1.0)
    grid = RectilinearGrid(size = N,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_capacity = relation.phase_transitions.density *
                    relation.phase_transitions.heat_capacity
    conductivity = diffusivity * heat_capacity
    energy_transport = ConductiveTemperatureTransport(conductivity = conductivity)
    heat_boundary_conditions = insulating_boundaries()

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        heat_boundary_conditions)

    return thermodynamics, heat_capacity
end

function initialize_cosine_energy_mode!(thermodynamics, heat_capacity;
                                        base_energy = -2e7,
                                        amplitude = 1e5)
    set!(thermodynamics;
         bulk_salinity = 0.0,
         temperature = z -> (base_energy + amplitude * cos(pi * z)) / heat_capacity)

    return nothing
end

function spatial_manufactured_solution_error(N; dt = 0.01, diffusivity = 1.0)
    thermodynamics, heat_capacity = pure_ice_conductive_column(N; diffusivity)
    base_energy = -2e7
    amplitude = 1e5

    initialize_cosine_energy_mode!(thermodynamics, heat_capacity;
                                   base_energy,
                                   amplitude)

    column_step!(thermodynamics, insulating_fluxes(), dt)

    numerical_energy = column_values(thermodynamics.fields.internal_energy)
    exact_amplitude = amplitude / (1 + diffusivity * pi^2 * dt)
    expected_energy = [base_energy + exact_amplitude * cos(pi * (k - 0.5) / N)
                       for k in 1:N]

    return maximum(abs.(numerical_energy .- expected_energy))
end

function temporal_manufactured_solution_error(dt; N = 128,
                                              final_time = 0.16,
                                              diffusivity = 1.0)
    thermodynamics, heat_capacity = pure_ice_conductive_column(N; diffusivity)
    base_energy = -2e7
    amplitude = 1e5
    steps = round(Int, final_time / dt)

    initialize_cosine_energy_mode!(thermodynamics, heat_capacity;
                                   base_energy,
                                   amplitude)

    for _ in 1:steps
        column_step!(thermodynamics, insulating_fluxes(), dt)
    end

    numerical_energy = column_values(thermodynamics.fields.internal_energy)
    discrete_eigenvalue = 4 * N^2 * sin(pi / (2N))^2
    exact_amplitude = amplitude * exp(-diffusivity * discrete_eigenvalue * final_time)
    expected_energy = [base_energy + exact_amplitude * cos(pi * (k - 0.5) / N)
                       for k in 1:N]

    return maximum(abs.(numerical_energy .- expected_energy))
end


salinity_variance(thermodynamics) = begin
    salinity = column_values(thermodynamics.fields.bulk_salinity)
    salinity_mean = sum(salinity) / length(salinity)
    sum(abs2, salinity .- salinity_mean) / length(salinity)
end

# Only the flux in excess of the start-of-step complete-melt energy ablates the surface, once per full step.
@testset "set! syncs the column grid metric" begin
    grid = RectilinearGrid(size = 8, z = SeaIceColumnDiscretization((0, 2)), topology = (Flat, Flat, Bounded))
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = IceWaterThermalEquilibrium(salinity = 0)))
    model = SeaIceModel(grid; ice_thermodynamics = thermodynamics)

    set!(model, h = 0.3, ℵ = 1)
    @test isapprox(column_height(1, 1, grid), 0.3; rtol = 1e-12)
end

