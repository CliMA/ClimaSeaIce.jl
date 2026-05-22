using ClimaSeaIce.SeaIceThermodynamics:
    BrineSalinityDiffusion,
    BulkSalinityDiffusion,
    ColumnBoundaryConditions,
    ColumnEnergyThermodynamics,
    ConductiveAndDiffusiveEnergyTransport,
    ConductiveTemperatureTransport,
    ExponentialShortwaveAbsorption,
    InsulatingBoundary,
    MeltingLimitedSurfaceFlux,
    NoSalinityTransport,
    NoShortwaveAbsorption,
    PrescribedBulkSalinity,
    PrescribedEnergyFlux,
    PrognosticBulkSalinity,
    QuadraticLiquidusEnergyRelation,
    assemble_column_energy_system!,
    brine_salinity,
    column_energy_budget,
    column_energy_time_step!,
    column_integrated_energy,
    column_integrated_salinity,
    column_salinity_time_step!,
    column_salt_budget,
    column_surface_stefan_residual_flux,
    column_stefan_thickness_budget,
    column_stefan_thickness_change,
    column_stefan_thickness_update!,
    compute_column_surface_stefan_residual_flux!,
    complete_melt_energy,
    compute_column_shortwave_flux!,
    compute_column_thermodynamic_diagnostics!,
    compute_column_transport_coefficients!,
    evolving_salinity_mushy_thermodynamics,
    internal_energy,
    liquid_fraction,
    prescribed_salinity_enthalpy_thermodynamics,
    solve_column_energy_system!,
    temperature,
    temperature_energy_derivative,
    temperature_salinity_derivative

using ClimaSeaIce: SeaIceModel
using Adapt
using Oceananigans
using Oceananigans.Fields: interior, set!
using Oceananigans: fields, prognostic_fields, prognostic_state, restore_prognostic_state!
using Statistics: median
using Test

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
    m = relation.phase_transitions.liquidus.slope
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
    m = liquidus.slope
    salinities, offsets = physical_column_relation_test_states(FT)

    max_error = zero(FT)
    for S in salinities, offset in offsets
        T = T0 - m * S - offset
        E = internal_energy(relation, T, S)
        brine_S = brine_salinity(relation, E, S)
        recovered_T = temperature(relation, E, S)
        liquidus_T = Float64(liquidus.freshwater_melting_temperature) -
                     Float64(liquidus.slope) * Float64(brine_S)
        max_error = max(max_error, abs(Float64(recovered_T) - liquidus_T))
    end

    return max_error
end

function test_liquid_fraction_bounds(::Type{FT}) where FT
    relation = QuadraticLiquidusEnergyRelation(FT)
    T0 = relation.phase_transitions.reference_temperature
    m = relation.phase_transitions.liquidus.slope
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
    m = relation.phase_transitions.liquidus.slope
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

@testset "Column energy thermodynamic relation" begin
    @testset "Float64 relation metrics" begin
        @test test_temperature_inversion(Float64) < 1e-11
        @test test_temperature_derivatives(Float64) < 1e-6

        temperature_error, energy_error, derivative_error = test_pure_ice_limit(Float64)
        @test temperature_error < 1e-12
        @test energy_error < 1e-8
        @test derivative_error < 1e-15

        @test test_liquidus_consistency(Float64) < 1e-12
        test_liquid_fraction_bounds(Float64)
        @test test_complete_melt_energy(Float64) < 1e-12
    end

    @testset "Float32 relation metrics" begin
        @test test_temperature_inversion(Float32) < 1e-5
        @test test_temperature_derivatives(Float32) < 1e-3

        temperature_error, energy_error, derivative_error = test_pure_ice_limit(Float32)
        @test temperature_error < 1e-6
        @test energy_error < 1
        @test derivative_error < 1e-7

        @test test_liquidus_consistency(Float32) < 1e-6
        test_liquid_fraction_bounds(Float32)
        @test test_complete_melt_energy(Float32) < 1e-6
    end
end

column_values(field) = vec(Array(interior(field)))

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
                           z = (0, 1),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_capacity = relation.phase_transitions.density *
                    relation.phase_transitions.heat_capacity
    conductivity = diffusivity * heat_capacity
    energy_transport = ConductiveTemperatureTransport(conductivity = conductivity)
    boundary_conditions = ColumnBoundaryConditions(top = InsulatingBoundary(),
                                                   bottom = InsulatingBoundary())

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        boundary_conditions)

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

    column_energy_time_step!(thermodynamics, dt)

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
        column_energy_time_step!(thermodynamics, dt)
    end

    numerical_energy = column_values(thermodynamics.fields.internal_energy)
    discrete_eigenvalue = 4 * N^2 * sin(pi / (2N))^2
    exact_amplitude = amplitude * exp(-diffusivity * discrete_eigenvalue * final_time)
    expected_energy = [base_energy + exact_amplitude * cos(pi * (k - 0.5) / N)
                       for k in 1:N]

    return maximum(abs.(numerical_energy .- expected_energy))
end

function performance_metric_column(size; prognostic_salinity = false)
    grid = RectilinearGrid(size = size,
                           z = (0, 1),
                           topology = (Flat, Flat, Bounded))

    boundary_conditions = ColumnBoundaryConditions(top = InsulatingBoundary(),
                                                   bottom = InsulatingBoundary())
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0)

    thermodynamics = if prognostic_salinity
        evolving_salinity_mushy_thermodynamics(grid;
            energy_transport,
            salinity_transport = BulkSalinityDiffusion(diffusivity = 1e-4),
            boundary_conditions)
    else
        prescribed_salinity_enthalpy_thermodynamics(grid;
            salinity_profile = 5.0,
            energy_transport,
            boundary_conditions)
    end

    set!(thermodynamics;
         bulk_salinity = z -> 5 + z,
         temperature = z -> -10 + z)

    for _ in 1:10
        column_energy_time_step!(thermodynamics, 1.0)
    end

    return thermodynamics
end

function warmed_column_step_allocations(; prognostic_salinity = false)
    thermodynamics = performance_metric_column(16; prognostic_salinity)
    allocations = [@allocated column_energy_time_step!(thermodynamics, 1.0)
                   for _ in 1:10]

    return maximum(allocations)
end

function median_column_step_runtime(size; steps = 30, samples = 5)
    thermodynamics = performance_metric_column(size)
    runtimes = zeros(samples)

    for n in 1:samples
        runtimes[n] = @elapsed begin
            for _ in 1:steps
                column_energy_time_step!(thermodynamics, 1.0)
            end
        end

        runtimes[n] /= steps
    end

    return median(runtimes)
end

@testset "Column energy thermodynamics containers" begin
    grid = RectilinearGrid(size = 4,
                           z = (0, 1),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    boundary_conditions = ColumnBoundaryConditions(top = PrescribedEnergyFlux(-20.0),
                                                   bottom = InsulatingBoundary())
    shortwave = ExponentialShortwaveAbsorption(surface_transmission = 0.3,
                                               attenuation_scale = 0.4)

    fixed = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 5.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
        shortwave_absorption = shortwave,
        boundary_conditions)

    @test fixed isa ColumnEnergyThermodynamics
    @test fixed.salinity_closure isa PrescribedBulkSalinity
    @test fixed.salinity_transport isa NoSalinityTransport
    @test fixed.shortwave_absorption === shortwave
    @test keys(fields(fixed)) == (:E,
                                  :bulk_salinity,
                                  :T,
                                  :liquid_fraction,
                                  :brine_salinity,
                                  :temperature_energy_derivative,
                                  :temperature_salinity_derivative)
    @test keys(prognostic_fields(fixed)) == (:E,)
    @test all(interior(fixed.fields.bulk_salinity) .== 5.0)

    set!(fixed; bulk_salinity = 4.0,
                temperature = -5.0)

    E = internal_energy(relation, -5.0, 4.0)
    phi = liquid_fraction(relation, E, 4.0)
    brine_S = brine_salinity(relation, E, 4.0)

    @test all(interior(fixed.fields.internal_energy) .≈ E)
    @test all(interior(fixed.fields.bulk_salinity) .== 4.0)
    @test all(interior(fixed.fields.temperature) .≈ -5.0)
    @test all(interior(fixed.fields.liquid_fraction) .≈ phi)
    @test all(interior(fixed.fields.brine_salinity) .≈ brine_S)

    evolving = evolving_salinity_mushy_thermodynamics(grid;
        relation,
        energy_transport = ConductiveAndDiffusiveEnergyTransport(conductivity = 2.0,
                                                                 diffusivity = 1e-9),
        salinity_transport = BulkSalinityDiffusion(diffusivity = 2e-9),
        shortwave_absorption = NoShortwaveAbsorption(),
        boundary_conditions)

    @test evolving.salinity_closure isa PrognosticBulkSalinity
    @test evolving.energy_transport isa ConductiveAndDiffusiveEnergyTransport
    @test evolving.salinity_transport isa BulkSalinityDiffusion
    @test keys(prognostic_fields(evolving)) == (:E, :bulk_salinity)
    @test Adapt.adapt(CPU(), evolving) isa ColumnEnergyThermodynamics

    set!(evolving; internal_energy = 3.0, bulk_salinity = 6.0)
    state = deepcopy(prognostic_state(evolving))
    set!(evolving; internal_energy = 4.0, bulk_salinity = 7.0)
    restore_prognostic_state!(evolving, state)

    @test all(interior(evolving.fields.internal_energy) .== 3.0)
    @test all(interior(evolving.fields.bulk_salinity) .== 6.0)

    model = SeaIceModel(grid;
                        ice_thermodynamics = evolving,
                        top_heat_flux = 0,
                        bottom_heat_flux = 0)

    @test model.ice_thermodynamics === evolving
    @test :E in keys(prognostic_fields(model))
    @test :bulk_salinity in keys(prognostic_fields(model))
    @test fields(model).E === evolving.fields.internal_energy
    @test fields(model).bulk_salinity === evolving.fields.bulk_salinity
end

@testset "Column energy fixed-grid time step" begin
    grid = RectilinearGrid(size = 8,
                           z = (0, 1),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    boundary_conditions = ColumnBoundaryConditions(top = InsulatingBoundary(),
                                                   bottom = InsulatingBoundary())
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0)
    dt = 5e3

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        boundary_conditions)

    set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> -12 + 4z)

    initial_energy = column_integrated_energy(thermodynamics)
    initial_profile = column_values(thermodynamics.fields.internal_energy)

    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_transport_coefficients!(thermodynamics)
    assemble_column_energy_system!(thermodynamics, dt)

    matrix = dense_column_energy_matrix(thermodynamics)
    rhs = column_values(thermodynamics.auxiliary.energy_rhs)
    expected_energy = matrix \ rhs

    solve_column_energy_system!(thermodynamics)
    stepped_energy = column_values(thermodynamics.fields.internal_energy)

    @test maximum(abs.(stepped_energy .- expected_energy)) < 1e-7
    @test abs(column_integrated_energy(thermodynamics) - initial_energy) < 1e-7
    @test maximum(abs.(stepped_energy .- initial_profile)) > 1e5

    constant = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        boundary_conditions)

    set!(constant; bulk_salinity = 0.0, temperature = -10.0)
    constant_energy = column_integrated_energy(constant)
    constant_profile = column_values(constant.fields.internal_energy)

    column_energy_time_step!(constant, dt)

    @test maximum(abs.(column_values(constant.fields.internal_energy) .- constant_profile)) < 1e-7
    @test abs(column_integrated_energy(constant) - constant_energy) < 1e-7
    @test maximum(abs.(column_values(constant.fields.temperature) .+ 10)) < 1e-12

    flux_boundary_conditions = ColumnBoundaryConditions(top = PrescribedEnergyFlux(2.0),
                                                        bottom = PrescribedEnergyFlux(0.5))
    flux_forced = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        boundary_conditions = flux_boundary_conditions)

    set!(flux_forced; bulk_salinity = 0.0, temperature = -10.0)
    flux_forced_energy = column_integrated_energy(flux_forced)
    column_energy_time_step!(flux_forced, 10.0)
    flux_budget = column_energy_budget(flux_forced, flux_forced_energy, 10.0)

    @test abs(column_integrated_energy(flux_forced) -
              flux_forced_energy -
              10.0 * (2.0 - 0.5)) < 1e-8
    @test abs(flux_budget.residual) < 1e-8
    @test flux_budget.relative_residual < 1e-11

    shortwave = ExponentialShortwaveAbsorption(surface_transmission = 3.0,
                                               attenuation_scale = 0.25)
    shortwave_forced = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        shortwave_absorption = shortwave,
        boundary_conditions)

    set!(shortwave_forced; bulk_salinity = 0.0, temperature = -10.0)
    shortwave_forced_energy = column_integrated_energy(shortwave_forced)
    compute_column_shortwave_flux!(shortwave_forced)

    @test column_values(shortwave_forced.auxiliary.shortwave_flux)[end] ≈ 3.0
    @test column_values(shortwave_forced.auxiliary.shortwave_flux)[1] ≈ 3.0 * exp(-4)

    column_energy_time_step!(shortwave_forced, 100.0)
    shortwave_budget = column_energy_budget(shortwave_forced, shortwave_forced_energy, 100.0)
    expected_shortwave_change = 100.0 * 3.0 * (1 - exp(-4))

    @test abs(shortwave_budget.shortwave_flux_change - expected_shortwave_change) < 1e-12
    @test abs(shortwave_budget.residual) < 1e-8
    @test shortwave_budget.relative_residual < 1e-11
end

@testset "Column energy manufactured convergence" begin
    spatial_errors = [spatial_manufactured_solution_error(N) for N in (16, 32, 64)]
    spatial_slopes = log2.(spatial_errors[1:2] ./ spatial_errors[2:3])

    @test all(spatial_errors .< [30, 8, 2])
    @test all(spatial_slopes .> 1.8)

    temporal_errors = [temporal_manufactured_solution_error(dt) for dt in (0.04, 0.02, 0.01)]
    temporal_slopes = log2.(temporal_errors[1:2] ./ temporal_errors[2:3])

    @test all(temporal_errors .< [6000, 3200, 1700])
    @test all(temporal_slopes .> 0.9)
end

salinity_variance(thermodynamics) = begin
    salinity = column_values(thermodynamics.fields.bulk_salinity)
    salinity_mean = sum(salinity) / length(salinity)
    sum(abs2, salinity .- salinity_mean) / length(salinity)
end

@testset "Column salinity fixed-grid time step" begin
    grid = RectilinearGrid(size = 8,
                           z = (0, 1),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    boundary_conditions = ColumnBoundaryConditions(top = InsulatingBoundary(),
                                                   bottom = InsulatingBoundary())
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0)
    salinity_transport = BulkSalinityDiffusion(diffusivity = 1e-4)
    dt = 100

    evolving = evolving_salinity_mushy_thermodynamics(grid;
        relation,
        energy_transport,
        salinity_transport,
        boundary_conditions)

    set!(evolving;
         bulk_salinity = z -> 5 + sin(2pi * z),
         temperature = -10.0)

    initial_salt = column_integrated_salinity(evolving)
    initial_salinity = column_values(evolving.fields.bulk_salinity)
    initial_variance = salinity_variance(evolving)

    column_salinity_time_step!(evolving, dt)
    salt_budget = column_salt_budget(evolving, initial_salt, dt)

    @test abs(column_integrated_salinity(evolving) - initial_salt) < 1e-12
    @test abs(salt_budget.residual) < 1e-12
    @test salt_budget.relative_residual < 1e-12
    @test salinity_variance(evolving) <= initial_variance + 1e-13
    @test maximum(abs.(column_values(evolving.fields.bulk_salinity) .- initial_salinity)) > 1e-3

    fixed = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        boundary_conditions)

    no_transport = evolving_salinity_mushy_thermodynamics(grid;
        relation,
        energy_transport,
        salinity_transport = NoSalinityTransport(),
        boundary_conditions)

    brine_marker = evolving_salinity_mushy_thermodynamics(grid;
        relation,
        energy_transport,
        salinity_transport = BrineSalinityDiffusion(diffusivity = 1e-4),
        boundary_conditions)

    set!(fixed;
         bulk_salinity = z -> 4 + z,
         temperature = z -> -15 + 3z)

    set!(no_transport;
         bulk_salinity = z -> 4 + z,
         temperature = z -> -15 + 3z)

    initial_no_transport_salinity = column_values(no_transport.fields.bulk_salinity)

    for _ in 1:100
        column_energy_time_step!(fixed, 1000)
        column_energy_time_step!(no_transport, 1000)
    end

    @test maximum(abs.(column_values(fixed.fields.internal_energy) .-
                       column_values(no_transport.fields.internal_energy))) < 1e-10
    @test maximum(abs.(column_values(fixed.fields.temperature) .-
                       column_values(no_transport.fields.temperature))) < 1e-10
    @test maximum(abs.(column_values(no_transport.fields.bulk_salinity) .-
                       initial_no_transport_salinity)) == 0
    @test maximum(abs.(column_values(fixed.fields.bulk_salinity) .-
                       initial_no_transport_salinity)) == 0
    @test_throws ArgumentError column_salinity_time_step!(brine_marker, dt)
end

@testset "Column Stefan thickness update" begin
    grid = RectilinearGrid(size = (3, 2, 1),
                           x = (0, 1),
                           y = (0, 1),
                           z = (0, 1),
                           topology = (Periodic, Periodic, Bounded))

    h = Field{Center, Center, Nothing}(grid)
    ρi = Field{Center, Center, Nothing}(grid)
    set!(h, 1.0)
    set!(ρi, 900.0)

    relation = QuadraticLiquidusEnergyRelation(Float64)
    phase_transitions = relation.phase_transitions
    residual_flux = 12.0
    dt = 3600.0
    expected_change = column_stefan_thickness_change(phase_transitions,
                                                     900.0,
                                                     residual_flux,
                                                     dt)

    column_stefan_thickness_update!(h,
                                    phase_transitions,
                                    ρi,
                                    residual_flux,
                                    dt)

    updated_thickness = Array(interior(h))
    budget = column_stefan_thickness_budget(1.0,
                                            1.0 + expected_change,
                                            phase_transitions,
                                            900.0,
                                            residual_flux,
                                            dt)

    @test maximum(abs.(updated_thickness .- (1.0 + expected_change))) < 1e-12
    @test abs(budget.residual) < 1e-15
    @test budget.relative_residual < 1e-12
end

@testset "Column melting-limited surface balance" begin
    grid = RectilinearGrid(size = 1,
                           z = (0, 1),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    bulk_salinity = 5.0
    initial_temperature = -2.0
    dt = 100.0
    excess_flux = 3.0
    initial_cell_energy = internal_energy(relation, initial_temperature, bulk_salinity)
    melt_cell_energy = complete_melt_energy(relation, bulk_salinity)
    flux_to_melt = (melt_cell_energy - initial_cell_energy) / dt
    requested_flux = flux_to_melt + excess_flux

    boundary_conditions = ColumnBoundaryConditions(top = MeltingLimitedSurfaceFlux(flux = requested_flux),
                                                   bottom = InsulatingBoundary())
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = bulk_salinity,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        boundary_conditions)

    set!(thermodynamics;
         bulk_salinity,
         temperature = initial_temperature)

    initial_energy = column_integrated_energy(thermodynamics)
    residual_flux = Field{Center, Center, Nothing}(grid)
    compute_column_surface_stefan_residual_flux!(residual_flux, thermodynamics, dt)
    scalar_residual_flux = column_surface_stefan_residual_flux(thermodynamics, dt)

    @test abs(first(interior(residual_flux)) + excess_flux) < 1e-12
    @test abs(scalar_residual_flux + excess_flux) < 1e-12

    column_energy_time_step!(thermodynamics, dt)
    budget = column_energy_budget(thermodynamics,
                                  initial_energy,
                                  dt;
                                  surface_stefan_residual_flux = scalar_residual_flux)

    @test abs(first(interior(thermodynamics.fields.internal_energy)) - melt_cell_energy) < 1e-7
    @test abs(column_integrated_energy(thermodynamics) - melt_cell_energy) < 1e-7
    @test abs(budget.surface_stefan_residual_change + dt * excess_flux) < 1e-12
    @test abs(budget.residual) < 1e-7
    @test budget.relative_residual < 1e-12

    h = Field{Center, Center, Nothing}(grid)
    set!(h, 1.0)
    ρi = 900.0
    expected_thickness = 1.0 + column_stefan_thickness_change(relation.phase_transitions,
                                                             ρi,
                                                             scalar_residual_flux,
                                                             dt)

    column_stefan_thickness_update!(h,
                                    relation.phase_transitions,
                                    ρi,
                                    residual_flux,
                                    dt)

    @test abs(first(interior(h)) - expected_thickness) < 1e-12
end

@testset "Column CPU performance metrics" begin
    @test warmed_column_step_allocations() < 2^11
    @test warmed_column_step_allocations(prognostic_salinity = true) < 2^11

    small_runtime = median_column_step_runtime(64; steps = 1000)
    large_runtime = median_column_step_runtime(128; steps = 1000)
    scaling_ratio = large_runtime / small_runtime

    @test 1.8 < scaling_ratio < 2.3
end
