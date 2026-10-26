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
    initialize_column_interfaces!,
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
using ClimaSeaIce: IncrementalRemapping, UniformThicknessDistribution, snow_slab_thermodynamics
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

@testset "Resolved ice and snow with a live atmospheric flux" begin
    grid = RectilinearGrid(Float64; size = (1, 1, 8), x = (0, 1), y = (0, 1),
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Bounded, Bounded, Bounded))
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    ice = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation, salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = FluxBoundary()))
    snow = snow_slab_thermodynamics(grid)
    atmosphere(i, j, grid, T, clock, fields) = 80 * (T + 30)
    top_flux = FluxFunction(atmosphere; top_temperature_dependent = true)
    model = SeaIceModel(grid; ice_thermodynamics = ice, snow_thermodynamics = snow,
                        phase_transitions = relation.phase_transitions,
                        top_heat_flux = top_flux, bottom_heat_flux = 0.0,
                        timestepper = :ForwardEuler)
    set!(model; h = 1.3, ℵ = 1.0, hs = 0.01, u = 0.0, v = 0.0)
    set!(ice; bulk_salinity = 0.0, temperature = -20.0)
    set!(snow.top_surface_temperature, -20.0)

    for _ in 1:12
        time_step!(model, 600.0)
        Ts = snow.top_surface_temperature[1, 1, 1]
        Ti = ice.auxiliary.surface_temperature[1, 1, 1]
        hs = model.snow_thickness[1, 1, 1]
        @test -30 < Ts < -20
        @test isapprox(0.31 * (Ti - Ts) / hs, atmosphere(1, 1, grid, Ts, model.clock, fields(model)); atol = 1e-5)
    end
end

@testset "Column properties follow conservative ice-volume remapping" begin
    N = 16
    grid = RectilinearGrid(Float64; size = (N, N, 3), x = (0, 1), y = (0, 1),
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Periodic, Periodic, Bounded), halo = (4, 4, 4))
    thermodynamics = evolving_salinity_mushy_thermodynamics(grid;
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        salinity_transport = NoSalinityTransport(),
        heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary()))
    model = SeaIceModel(grid; ice_thermodynamics = thermodynamics,
                        advection = IncrementalRemapping(), timestepper = :ForwardEuler,
                        top_heat_flux = 0.0, bottom_heat_flux = 0.0)
    anomaly(x, y) = exp(-((x - 0.35)^2 + (y - 0.45)^2) / 0.1^2)
    set!(model; h = 1.0, ℵ = (x, y) -> 0.5anomaly(x, y),
         u = (x, y) -> 0.2, v = (x, y) -> 0.1)
    set!(thermodynamics; bulk_salinity = (x, y, z) -> 4 + anomaly(x, y),
         temperature = (x, y, z) -> -5 - anomaly(x, y))

    E⁰ = Array(interior(thermodynamics.fields.internal_energy))
    V⁰ = Array(interior(model.ice_thickness)) .* Array(interior(model.ice_concentration))
    S⁰ = Array(interior(thermodynamics.fields.bulk_salinity))
    energy⁰ = sum(E⁰ .* V⁰)
    salt⁰ = sum(S⁰ .* V⁰)

    for _ in 1:20
        time_step!(model, 0.05)
    end

    E = Array(interior(thermodynamics.fields.internal_energy))
    S = Array(interior(thermodynamics.fields.bulk_salinity))
    V = Array(interior(model.ice_thickness)) .* Array(interior(model.ice_concentration))
    @test isapprox(sum(V), sum(V⁰); rtol = 1e-12)
    @test isapprox(sum(E .* V), energy⁰; rtol = 1e-12)
    @test isapprox(sum(S .* V), salt⁰; rtol = 1e-12)
    @test maximum(abs.(E .- E⁰)) > 1e-6
    @test all(isfinite, E) && all(isfinite, S)
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
    m = relation.phase_transitions.liquidus.slope
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
    m = phase_transitions.liquidus.slope
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

@testset "Fixed drained-ice salinity profile" begin
    profile = FixedDrainedIceSalinityProfile(Float64)
    depths = [(2k - 1) // 14 for k in 1:7]
    expected = Float64[bl99_salinity_reference(s) for s in depths]
    actual = [salinity_at_normalized_depth(profile, Float64(s)) for s in depths]

    @test maximum(abs.(actual .- expected)) < 1e-14
    @test abs(salinity_at_normalized_depth(profile, 0.0)) < 1e-14
    @test abs(salinity_at_normalized_depth(profile, 1.0) - profile.maximum_salinity) < 1e-14

    zetas = range(0, 1, length = 17)
    reversal_errors = [abs(profile(zeta) -
                       salinity_at_normalized_depth(profile, 1 - zeta))
                       for zeta in zetas]
    height_errors = [abs(salinity_at_normalized_height(profile, zeta) -
                         salinity_at_normalized_depth(profile, 1 - zeta))
                     for zeta in zetas]

    @test maximum(reversal_errors) < 1e-14
    @test maximum(height_errors) < 1e-14

    # A regional or global ice column has two active horizontal coordinates.
    grid = RectilinearGrid(Float64; size = (2, 2, 8), x = (0, 1), y = (0, 1),
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Bounded, Bounded, Bounded))
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid; salinity_profile = profile)
    salinity = interior(thermodynamics.fields.bulk_salinity)
    for k in 1:8
        expected = profile((k - 0.5) / 8)
        @test maximum(abs.(salinity[:, :, k] .- expected)) < 1e-14
    end
    @test salinity[1, 1, 1] > salinity[1, 1, 8] # saltier at the ice base

    set!(thermodynamics; bulk_salinity = 0.0)
    set!(thermodynamics; bulk_salinity = profile)
    salinity = interior(thermodynamics.fields.bulk_salinity)
    @test maximum(abs.(salinity[:, :, 8] .- profile(7.5 / 8))) < 1e-14

    # The profile is tied to normalized ice depth, even after the physical
    # column interfaces move from [0, 1] to [-2, 0].
    ice_thickness = Field{Center, Center, Nothing}(grid)
    set!(ice_thickness, 2.0)
    initialize_column_interfaces!(grid, ice_thickness)
    set!(thermodynamics; bulk_salinity = profile)
    salinity = interior(thermodynamics.fields.bulk_salinity)
    for k in 1:8
        @test maximum(abs.(salinity[:, :, k] .- profile((k - 0.5) / 8))) < 1e-14
    end

    # Grid active-cell maps skip masked cells during :xyz launches. Initial ice
    # data can still occupy those cells, so the profile setter fills the field.
    inactive_first_cell(x, y, z) = x < 0.5 && y < 0.5
    masked_grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(inactive_first_cell); active_cells_map = true)
    masked_thermodynamics = prescribed_salinity_enthalpy_thermodynamics(masked_grid; salinity_profile = profile)
    masked_salinity = interior(masked_thermodynamics.fields.bulk_salinity)
    for k in 1:8
        @test maximum(abs.(masked_salinity[:, :, k] .- profile((k - 0.5) / 8))) < 1e-14
    end
end

@testset "Fixed salinity brine-pocket relation" begin
    @testset "Float64 relation metrics" begin
        @test test_fixed_salinity_temperature_inversion(Float64) < 1e-11
        @test test_fixed_salinity_temperature_derivatives(Float64) < 1e-6

        temperature_error, energy_error, derivative_error =
            test_fixed_salinity_pure_ice_limit(Float64)
        @test temperature_error < 1e-12
        @test energy_error < 1e-8
        @test derivative_error < 1e-15

        @test test_fixed_salinity_complete_melt_threshold(Float64) < 1e-11
        @test test_fixed_salinity_heat_capacity_identity(Float64) < 1e-12
    end

    @testset "Float32 relation metrics" begin
        @test test_fixed_salinity_temperature_inversion(Float32) < 1e-5
        @test test_fixed_salinity_temperature_derivatives(Float32) < 1e-3

        temperature_error, energy_error, derivative_error =
            test_fixed_salinity_pure_ice_limit(Float32)
        @test temperature_error < 1e-6
        @test energy_error < 1
        @test derivative_error < 1e-7

        @test test_fixed_salinity_complete_melt_threshold(Float32) < 1e-5
        @test test_fixed_salinity_heat_capacity_identity(Float32) < 1e-6
    end
end

@testset "CICE-compatible thermal conductivity closures" begin
    maykut = MaykutUntersteinerConductivity(Float64)
    bubbly = BubblyBrineConductivity(Float64)
    states = ((-20.0, 0.0), (-5.0, 3.2), (-2.0, 5.0), (-1e-12, 5.0))

    for (T, S) in states
        @test isapprox(ice_thermal_conductivity(maykut, T, S),
                       maykut_untersteiner_reference(maykut, T, S);
                       atol = 1e-14, rtol = 0)
        @test isapprox(ice_thermal_conductivity(bubbly, T, S),
                       bubbly_brine_reference(bubbly, T, S);
                       atol = 1e-14, rtol = 0)
        @test isfinite(ice_thermal_conductivity(maykut, T, S))
        @test isfinite(ice_thermal_conductivity(bubbly, T, S))
        @test ice_thermal_conductivity(maykut, T, S) >= maykut.minimum_conductivity
        @test ice_thermal_conductivity(bubbly, T, S) >= bubbly.minimum_conductivity
    end

    lower_T, lower_S, lower_h = -10.0, 1.0, 0.1
    upper_T, upper_S, upper_h = -2.0, 5.0, 0.3

    lower_K = ice_thermal_conductivity(maykut, lower_T, lower_S)
    upper_K = ice_thermal_conductivity(maykut, upper_T, upper_S)
    expected_face_K = (lower_h + upper_h) / (lower_h / lower_K + upper_h / upper_K)

    @test isapprox(face_thermal_conductivity(maykut,
                                             lower_T, lower_S, lower_h,
                                             upper_T, upper_S, upper_h),
                   expected_face_K;
                   atol = 1e-14, rtol = 0)

    @test face_thermal_conductivity(2.0,
                                    lower_T, lower_S, lower_h,
                                    upper_T, upper_S, upper_h) === 2.0

    grid = RectilinearGrid(size = 4,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    profile = FixedDrainedIceSalinityProfile(Float64)
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    energy_transport = ConductiveTemperatureTransport(conductivity = maykut)
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = profile,
        energy_transport,
        heat_boundary_conditions = insulating_boundaries())

    set!(thermodynamics;
         bulk_salinity = profile,
         temperature = z -> -15 + 10z)

    compute_column_transport_coefficients!(thermodynamics)

    temperatures = column_values(thermodynamics.fields.temperature)
    salinities = column_values(thermodynamics.fields.bulk_salinity)
    face_conductivities = column_values(thermodynamics.auxiliary.thermal_conductivity)
    dz = 0.25

    @test isapprox(face_conductivities[1],
                   ice_thermal_conductivity(maykut, temperatures[1], salinities[1]);
                   atol = 1e-14, rtol = 0)
    @test isapprox(face_conductivities[end],
                   ice_thermal_conductivity(maykut, temperatures[end], salinities[end]);
                   atol = 1e-14, rtol = 0)

    for k in 2:4
        expected_K = face_thermal_conductivity(maykut,
                                               temperatures[k-1],
                                               salinities[k-1],
                                               dz,
                                               temperatures[k],
                                               salinities[k],
                                               dz)
        @test isapprox(face_conductivities[k], expected_K; atol = 1e-14, rtol = 0)
    end
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


@testset "Column energy thermodynamics containers" begin
    grid = RectilinearGrid(size = 4,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary())
    shortwave = ExponentialShortwaveAbsorption(surface_transmission = 0.3,
                                               attenuation_scale = 0.4)

    fixed = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 5.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
        shortwave_absorption = shortwave,
        heat_boundary_conditions)

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
        energy_transport = ConductiveTemperatureTransport(conductivity = 2.0, diffusivity = 1e-9),
        salinity_transport = BulkSalinityDiffusion(diffusivity = 2e-9),
        shortwave_absorption = NoShortwaveAbsorption(),
        heat_boundary_conditions)

    @test evolving.salinity_closure isa PrognosticBulkSalinity
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
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_boundary_conditions = insulating_boundaries()
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0)
    dt = 5e3

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        heat_boundary_conditions)

    set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> -12 + 4z)

    initial_energy = column_integrated_energy(thermodynamics)
    initial_profile = column_values(thermodynamics.fields.internal_energy)

    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_transport_coefficients!(thermodynamics)
    assemble_column_energy_system!(thermodynamics, insulating_fluxes(), TEST_CLOCK, NO_FIELDS, CONSOLIDATED, dt)

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
        heat_boundary_conditions)

    set!(constant; bulk_salinity = 0.0, temperature = -10.0)
    constant_energy = column_integrated_energy(constant)
    constant_profile = column_values(constant.fields.internal_energy)

    column_step!(constant, insulating_fluxes(), dt)

    @test maximum(abs.(column_values(constant.fields.internal_energy) .- constant_profile)) < 1e-7
    @test abs(column_integrated_energy(constant) - constant_energy) < 1e-7
    @test maximum(abs.(column_values(constant.fields.temperature) .+ 10)) < 1e-12

    flux_fluxes = (top = 2.0, bottom = 0.5)
    flux_forced = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary()))

    set!(flux_forced; bulk_salinity = 0.0, temperature = -10.0)
    flux_forced_energy = column_integrated_energy(flux_forced)
    column_step!(flux_forced, flux_fluxes, 10.0)

    # Upward-positive fluxes: net energy gain is Δt·(Qᵇ − Qᵘ) = 10·(0.5 − 2.0).
    @test abs(column_integrated_energy(flux_forced) -
              flux_forced_energy -
              10.0 * (0.5 - 2.0)) < 1e-8

    shortwave = ExponentialShortwaveAbsorption(surface_transmission = 3.0,
                                               attenuation_scale = 0.25)
    shortwave_forced = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        shortwave_absorption = shortwave,
        heat_boundary_conditions)

    set!(shortwave_forced; bulk_salinity = 0.0, temperature = -10.0)
    shortwave_forced_energy = column_integrated_energy(shortwave_forced)
    compute_column_shortwave_flux!(shortwave_forced)

    @test column_values(shortwave_forced.auxiliary.shortwave_flux)[end] ≈ 3.0
    @test column_values(shortwave_forced.auxiliary.shortwave_flux)[1] ≈ 3.0 * exp(-4)

    column_step!(shortwave_forced, insulating_fluxes(), 100.0)
    expected_shortwave_change = 100.0 * 3.0 * (1 - exp(-4))

    @test abs(column_integrated_energy(shortwave_forced) - shortwave_forced_energy - expected_shortwave_change) < 1e-8

    bottom_temperature = -1.8
    top_flux = -3.0
    prescribed_bottom = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        heat_boundary_conditions = (top = FluxBoundary(),
                                    bottom = PrescribedTemperature(bottom_temperature)))

    set!(prescribed_bottom; bulk_salinity = 0.0, temperature = z -> -12 + 4z)

    initial_temperature = column_values(prescribed_bottom.fields.temperature)
    n = length(initial_temperature)
    Δz = 1 / n
    ρc = relation.phase_transitions.density * relation.phase_transitions.heat_capacity
    η = dt / (ρc * Δz)
    conductance = 2 * 2.0 / Δz
    interface_conductance = 2.0 / Δz
    lower_diagonal = zeros(n - 1)
    diagonal = ones(n)
    upper_diagonal = zeros(n - 1)
    rhs = copy(initial_temperature)

    for k in 1:n
        if k == 1
            diagonal[k] += η * (conductance + interface_conductance)
            upper_diagonal[k] = -η * interface_conductance
            rhs[k] += η * conductance * bottom_temperature
        elseif k == n
            lower_diagonal[k-1] = -η * interface_conductance
            diagonal[k] += η * interface_conductance
            rhs[k] -= η * top_flux  # upward-positive top flux enters the top cell as -Qᵘ
        else
            lower_diagonal[k-1] = -η * interface_conductance
            upper_diagonal[k] = -η * interface_conductance
            diagonal[k] += 2 * η * interface_conductance
        end
    end

    expected_temperature =
        solve_reference_vector_tridiagonal(lower_diagonal, diagonal, upper_diagonal, rhs)

    column_step!(prescribed_bottom, (top = top_flux, bottom = 0.0), dt)

    @test maximum(abs.(column_values(prescribed_bottom.fields.temperature) .-
                       expected_temperature)) < 1e-11
end

@testset "Column energy moving-grid metric" begin
    grid = RectilinearGrid(size = 4,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_boundary_conditions = insulating_boundaries()
    no_transport = ConductiveTemperatureTransport(conductivity = 0.0)

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = no_transport,
        heat_boundary_conditions)

    set!(thermodynamics; bulk_salinity = 0.0, temperature = -10.0)

    initial_energy = column_integrated_energy(thermodynamics)
    initial_profile = column_values(thermodynamics.fields.internal_energy)
    current_metric = 0.75

    set_column_metric!(grid, 1.0, current_metric)
    column_step!(thermodynamics, insulating_fluxes(), 10.0)

    @test column_integrated_energy(thermodynamics) ≈ current_metric * initial_energy
    @test column_values(thermodynamics.fields.internal_energy) ≈ initial_profile

    set_column_metric!(grid, 1.0, 1.0)
    set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> -12 + 8z)
    nonuniform_initial_profile = column_values(thermodynamics.fields.internal_energy)
    source_faces = collect(range(0.0, 1.0; length = 5))
    contracted_faces = collect(range(0.0, 0.75; length = 5))
    expected_contracted_profile =
        conservative_column_remap(nonuniform_initial_profile, source_faces, contracted_faces)

    set_column_metric!(grid, 1.0, 0.75)
    column_step!(thermodynamics, insulating_fluxes(), 10.0)

    @test column_values(thermodynamics.fields.internal_energy) ≈ expected_contracted_profile

    set_column_metric!(grid, 1.0, 1.0)
    set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> -12 + 8z)
    nonuniform_initial_profile = column_values(thermodynamics.fields.internal_energy)
    expanded_faces = collect(range(0.0, 1.25; length = 5))
    expected_expanded_profile =
        conservative_column_remap(nonuniform_initial_profile,
                                  source_faces,
                                  expanded_faces;
                                  fill_value = last(nonuniform_initial_profile))

    set_column_metric!(grid, 1.0, 1.25)
    column_step!(thermodynamics, insulating_fluxes(), 10.0)

    @test column_values(thermodynamics.fields.internal_energy) ≈ expected_expanded_profile

    basal_grid = RectilinearGrid(size = 4,
                                 z = SeaIceColumnDiscretization((0, 1)),
                                 topology = (Flat, Flat, Bounded))
    # Congelation lays down ice at the basal Dirichlet temperature; the swept enthalpy is the BL99 new-ice
    # enthalpy internal_energy(relation, Tᵇ, S). With zero conductivity the Dirichlet face adds nothing to the
    # interior solve, isolating the conservative basal-fill remap.
    basal_bottom_temperature = -1.8
    fill_energy = internal_energy(relation, basal_bottom_temperature, 0.0)
    basal_growth = prescribed_salinity_enthalpy_thermodynamics(basal_grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = no_transport,
        heat_boundary_conditions = (top = FluxBoundary(),
                                    bottom = PrescribedTemperature(basal_bottom_temperature)))

    set_column_metric!(basal_grid, 1.0, 1.0; surface = 1.0, move = :base)
    set!(basal_growth; bulk_salinity = 0.0, temperature = z -> -12 + 8z)
    basal_initial_profile = column_values(basal_growth.fields.internal_energy)
    basal_growth_faces = collect(range(-0.25, 1.0; length = 5))
    expected_basal_growth_profile =
        conservative_column_remap(basal_initial_profile,
                                  source_faces,
                                  basal_growth_faces;
                                  fill_value = fill_energy)

    set_column_metric!(basal_grid, 1.0, 1.25; surface = 1.0, move = :base)
    column_step!(basal_growth, insulating_fluxes(), 10.0)

    @test column_values(basal_growth.fields.internal_energy) ≈ expected_basal_growth_profile

    flux_forced = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = no_transport,
        heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary()))

    set_column_metric!(grid, 1.0, 1.0)
    set!(flux_forced; bulk_salinity = 0.0, temperature = -10.0)
    flux_initial_energy = column_integrated_energy(flux_forced)

    set_column_metric!(grid, 1.0, 1.25)
    column_step!(flux_forced, (top = 2.0, bottom = 0.0), 10.0)

    # Upward-positive top flux removes energy: net change is Δt·(−Qᵘ) = −20.
    @test abs(column_integrated_energy(flux_forced) - 1.25 * flux_initial_energy + 20.0) < 1e-8

    salinity_grid = RectilinearGrid(size = 4,
                                    z = SeaIceColumnDiscretization((0, 1)),
                                    topology = (Flat, Flat, Bounded))

    evolving = evolving_salinity_mushy_thermodynamics(salinity_grid;
        relation,
        energy_transport = no_transport,
        salinity_transport = NoSalinityTransport(),
        heat_boundary_conditions)

    set!(evolving; bulk_salinity = 5.0, temperature = -10.0)
    initial_salt = column_integrated_salinity(evolving)
    initial_salinity = column_values(evolving.fields.bulk_salinity)

    set_column_metric!(salinity_grid, 1.0, 1.25)
    column_salinity_time_step!(evolving, 10.0)

    @test column_integrated_salinity(evolving) ≈ 1.25 * initial_salt
    @test column_values(evolving.fields.bulk_salinity) ≈ initial_salinity

    set_column_metric!(salinity_grid, 1.0, 1.0)
    set!(evolving; bulk_salinity = z -> 2 + 6z, temperature = -10.0)
    initial_salinity = column_values(evolving.fields.bulk_salinity)
    expected_salinity =
        conservative_column_remap(initial_salinity, source_faces, contracted_faces)

    set_column_metric!(salinity_grid, 1.0, 0.75)
    column_salinity_time_step!(evolving, 10.0)

    @test column_values(evolving.fields.bulk_salinity) ≈ expected_salinity
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
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_boundary_conditions = insulating_boundaries()
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0)
    salinity_transport = BulkSalinityDiffusion(diffusivity = 1e-4)
    dt = 100

    evolving = evolving_salinity_mushy_thermodynamics(grid;
        relation,
        energy_transport,
        salinity_transport,
        heat_boundary_conditions)

    set!(evolving;
         bulk_salinity = z -> 5 + sin(2pi * z),
         temperature = -10.0)

    initial_salt = column_integrated_salinity(evolving)
    initial_salinity = column_values(evolving.fields.bulk_salinity)
    initial_variance = salinity_variance(evolving)

    column_salinity_time_step!(evolving, dt)

    @test abs(column_integrated_salinity(evolving) - initial_salt) < 1e-12
    @test salinity_variance(evolving) <= initial_variance + 1e-13
    @test maximum(abs.(column_values(evolving.fields.bulk_salinity) .- initial_salinity)) > 1e-3

    fixed = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        heat_boundary_conditions)

    no_transport = evolving_salinity_mushy_thermodynamics(grid;
        relation,
        energy_transport,
        salinity_transport = NoSalinityTransport(),
        heat_boundary_conditions)

    set!(fixed;
         bulk_salinity = z -> 4 + z,
         temperature = z -> -15 + 3z)

    set!(no_transport;
         bulk_salinity = z -> 4 + z,
         temperature = z -> -15 + 3z)

    initial_no_transport_salinity = column_values(no_transport.fields.bulk_salinity)

    for _ in 1:100
        column_step!(fixed, insulating_fluxes(), 1000)
        column_step!(no_transport, insulating_fluxes(), 1000)
    end

    @test maximum(abs.(column_values(fixed.fields.internal_energy) .-
                       column_values(no_transport.fields.internal_energy))) < 1e-10
    @test maximum(abs.(column_values(fixed.fields.temperature) .-
                       column_values(no_transport.fields.temperature))) < 1e-10
    @test maximum(abs.(column_values(no_transport.fields.bulk_salinity) .-
                       initial_no_transport_salinity)) == 0
    @test maximum(abs.(column_values(fixed.fields.bulk_salinity) .-
                       initial_no_transport_salinity)) == 0
end

@testset "Melting surface partitions heat between conduction and ablation" begin
    grid = RectilinearGrid(size = 1, z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    S = 4.0
    T₀ = -5.0
    Δt = 100.0
    K = 2.0
    Qᵃ = -100.0 # upward positive, so this warms the ice
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation, salinity_profile = S,
        energy_transport = ConductiveTemperatureTransport(conductivity = K),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = FluxBoundary()))
    set!(thermodynamics; temperature = T₀)
    column_step!(thermodynamics, (top = Qᵃ, bottom = 0.0), Δt)

    Tₛ = first(interior(thermodynamics.auxiliary.surface_temperature))
    Tᵢ = first(interior(thermodynamics.fields.temperature))
    Tₘ = -relation.phase_transitions.liquidus.slope * S
    Qᶜ = 2K * (Tᵢ - Tₛ) # top-face conduction at unit column thickness
    Qˢ = column_surface_stefan_residual_flux(thermodynamics.heat_boundary_conditions.top,
                                             Qᵃ, 1, 1, 1, grid, thermodynamics.auxiliary,
                                             thermodynamics.fields, relation, TEST_CLOCK, NO_FIELDS, Δt)
    @test Tₛ ≈ Tₘ
    @test Qˢ < 0
    @test isapprox(Qˢ, Qᵃ - Qᶜ; atol = 1e-8)
    @test Tᵢ > T₀ # conductive heat enters the column while the surface melts
end

@testset "Surface ablation uses the enthalpy of the disappearing ice" begin
    grid = RectilinearGrid(size = 8, z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    S, T₀, Δt, Qᵃ = 4.0, -5.0, 900.0, -100.0
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation, salinity_profile = S,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = FluxBoundary()))
    model = SeaIceModel(grid; ice_thermodynamics = thermodynamics,
                        phase_transitions = relation.phase_transitions,
                        ice_consolidation_thickness = 0.05,
                        top_heat_flux = Qᵃ, bottom_heat_flux = 0.0)
    set!(model; h = 1.5, ℵ = 1.0)
    set!(thermodynamics; temperature = T₀)

    V⁰ = first(interior(model.ice_thickness)) * first(interior(model.ice_concentration))
    E⁰ = internal_energy(relation, T₀, S)
    Tₘ = -relation.phase_transitions.liquidus.slope * S
    𝒫 = relation.phase_transitions
    # Icepack's BL99 surface-layer enthalpy q_i(T, S), expressed as a positive
    # energy per volume required to melt ice and warm its meltwater to 0°C.
    melt_energy = 𝒫.density * (𝒫.heat_capacity * (Tₘ - T₀) +
                               𝒫.reference_latent_heat * (1 - Tₘ / T₀) -
                               𝒫.liquid_heat_capacity * Tₘ)
    time_step!(model, Δt)
    V¹ = first(interior(model.ice_thickness)) * first(interior(model.ice_concentration))
    E¹ = last(interior(thermodynamics.fields.internal_energy))

    @test isapprox(V¹ - V⁰, Qᵃ * Δt / melt_energy; rtol = 1e-6)
    @test isapprox(E¹, E⁰; rtol = 1e-12) # no conduction into the top cell
    @test E¹ < complete_melt_energy(relation, S) # ablation starts before the cell fully melts
    # The shifted ice energy includes its latent reference; use concentration
    # times thickness because proportional evolution changes both fields.
    ρL = 𝒫.density * 𝒫.reference_latent_heat
    @test isapprox(V¹ * (E¹ - ρL) - V⁰ * (E⁰ - ρL), -Qᵃ * Δt; atol = 1e-5)
end

@testset "set! syncs the column grid metric" begin
    grid = RectilinearGrid(size = 8, z = SeaIceColumnDiscretization((0, 2)), topology = (Flat, Flat, Bounded))
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = IceWaterThermalEquilibrium(salinity = 0)))
    model = SeaIceModel(grid; ice_thermodynamics = thermodynamics)

    set!(model, h = 0.3, ℵ = 1)
    @test isapprox(column_height(grid, 1, 1), 0.3; rtol = 1e-12)
end

@testset "Column implicit surface temperature solve" begin
    Nz = 8
    grid = RectilinearGrid(size = Nz,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    conductivity = 2.0
    base_temperature = -8.0
    dt = 1e3

    # Linearized warming atmospheric flux Q_atm(Tₛ) = Q₀ - λ (Tₛ - T_ref), evaluated at the surface temperature.
    # A T-dependent flux is what makes the outer iteration non-trivial — a scalar flux would converge in one pass.
    Q₀, λ, T_ref = 5.0, 0.5, -10.0
    atmospheric_flux(i, j, grid, T, clock, fields) = Q₀ - λ * (T - T_ref)
    top_flux = FluxFunction(atmospheric_flux; top_temperature_dependent = true)
    ext = (top = top_flux, bottom = 0.0)

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = conductivity),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(),
                                    bottom = PrescribedTemperature(base_temperature)))

    set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> base_temperature + 2z)

    for _ in 1:50
        column_step!(thermodynamics, ext, dt)
    end

    surface_temperature = first(interior(thermodynamics.auxiliary.surface_temperature))
    top_temperature = column_values(thermodynamics.fields.temperature)[end]
    Δz = 1 / Nz
    conductance = 2 * conductivity / Δz
    melt_temperature = 0.0

    @test surface_temperature <= melt_temperature + 1e-12
    # Upward-positive massless-surface balance: Q_atm(Tₛ) = conductance·(T_top − Tₛ), to the iteration tolerance.
    @test abs((Q₀ - λ * (surface_temperature - T_ref)) -
              conductance * (top_temperature - surface_temperature)) < 1e-6
end

@testset "Column with a prescribed top temperature builds without top_heat_flux" begin
    grid = RectilinearGrid(size = 4, z = SeaIceColumnDiscretization((0, 1)), topology = (Flat, Flat, Bounded))
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        heat_boundary_conditions = (top = PrescribedTemperature(-10.0), bottom = IceWaterThermalEquilibrium(salinity = 0)))
    model = SeaIceModel(grid; ice_thermodynamics = thermodynamics)

    set!(model, h = 0.5, ℵ = 1)
    set!(thermodynamics; bulk_salinity = 0.0, temperature = -5.0)
    time_step!(model, 600)

    @test model.external_heat_fluxes.top == 0
    @test first(interior(model.ice_thickness)) > 0.5
end

@testset "Resolved ice column with finite slab snow closes its surface budget" begin
    function snow_case(Qatm, hs⁰; Δt = 600.0, initial_temperature = -5.0,
                       h⁰ = 1.0, snowfall = 0.0, snow_itd_shape = nothing,
                       ice_conductivity = 2.0)
        Nz = 8
        grid = RectilinearGrid(Float64; size = (1, 1, Nz), x = (0, 1), y = (0, 1),
                               z = SeaIceColumnDiscretization((0, 1)),
                               topology = (Bounded, Bounded, Bounded))
        relation = FixedSalinityBrinePocketEnergyRelation(Float64)
        ice = prescribed_salinity_enthalpy_thermodynamics(grid; relation,
            salinity_profile = 4.0,
            energy_transport = ConductiveTemperatureTransport(conductivity = ice_conductivity),
            heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(),
                                        bottom = FluxBoundary()))
        snow = snow_slab_thermodynamics(grid; itd_shape = snow_itd_shape)
        model = SeaIceModel(grid; ice_thermodynamics = ice, snow_thermodynamics = snow,
                            phase_transitions = relation.phase_transitions,
                            top_heat_flux = Qatm, bottom_heat_flux = 0.0,
                            timestepper = :ForwardEuler, snowfall)
        set!(model; h = h⁰, ℵ = 1.0, hs = hs⁰, u = 0.0, v = 0.0)
        set!(ice; bulk_salinity = 4.0, temperature = initial_temperature)
        E⁰ = sum(interior(ice.fields.internal_energy)) * h⁰ / Nz
        time_step!(model, Δt)
        h = first(interior(model.ice_thickness))
        ℵ = first(interior(model.ice_concentration))
        hs = first(interior(model.snow_thickness))
        E = sum(interior(ice.fields.internal_energy)) * h / Nz
        Ts = first(interior(snow.top_surface_temperature))
        Ti = first(interior(ice.auxiliary.surface_temperature))
        ρs = model.snow_density[1, 1, 1]
        snow_melt_energy = ρs * model.phase_transitions.reference_latent_heat * (hs⁰ - hs)
        residual = E - E⁰ + snow_melt_energy + Qatm * Δt
        skin_residual = first(interior(ice.auxiliary.snow_surface_residual))
        ice_mass_flux = model.mass_fluxes.thermodynamics.ice[1, 1, 1]
        snow_mass_flux = model.mass_fluxes.thermodynamics.snow[1, 1, 1]
        intercepted_snowfall = model.mass_fluxes.intercepted_snowfall[1, 1, 1]
        return (; h, ℵ, hs, Ts, Ti, residual, skin_residual,
                ice_mass_flux, snow_mass_flux, intercepted_snowfall)
    end

    cold = snow_case(20.0, 0.2)
    distributed = snow_case(20.0, 0.2; snow_itd_shape = UniformThicknessDistribution(Float64; categories = 4))
    thin = snow_case(20.0, 1e-7)
    dust = snow_case(20.0, 1e-12)
    warm = snow_case(-100.0, 0.2)
    exhausted = snow_case(-2000.0, 1e-4; Δt = 20.0, initial_temperature = -0.3)
    bare_melt = snow_case(-100.0, 1e-12; Δt = 900.0, h⁰ = 1.5, ice_conductivity = 0.0)
    accumulating = snow_case(0.0, 0.1; snowfall = 0.003)
    flooded = snow_case(0.0, 0.2; h⁰ = 0.1)
    @test cold.hs == 0.2 && cold.Ts < 0
    @test distributed.Ts > cold.Ts
    @test isapprox(0.31 * UniformThicknessDistribution(Float64; categories = 4).conductivity_factor *
                   (distributed.Ti - distributed.Ts) / distributed.hs, 20.0; atol = 1e-5)
    @test distributed.skin_residual < 1e-6
    @test thin.hs ≈ 1e-7 && thin.Ts < 0 && thin.skin_residual < 1e-5
    @test dust.hs ≈ 1e-12 && dust.Ts == dust.Ti && dust.skin_residual == 0
    @test warm.hs < 0.2 && warm.Ts == 0
    @test exhausted.hs < 1e-7 && exhausted.Ts == 0 && exhausted.ice_mass_flux < 0
    phase_transitions = FixedSalinityBrinePocketEnergyRelation(Float64).phase_transitions
    Tₘ = -phase_transitions.liquidus.slope * 4.0
    melt_energy = phase_transitions.density * (
        phase_transitions.heat_capacity * (Tₘ + 5.0) +
        phase_transitions.reference_latent_heat * (1 + Tₘ / 5.0) -
        phase_transitions.liquid_heat_capacity * Tₘ)
    @test isapprox(bare_melt.h * bare_melt.ℵ - 1.5, -100 * 900 / melt_energy; rtol = 1e-6)
    # `residual` omits the enthalpy carried away by ablated ice, so it only
    # measures a closed ice-plus-snow budget in cases without ice ablation.
    @test maximum(abs.((cold.residual, warm.residual))) < 1e-6
    @test maximum((cold.skin_residual, warm.skin_residual, exhausted.skin_residual)) < 1e-6
    @test isapprox(accumulating.hs, 0.1 + 0.003 * 600 / 330; rtol = 1e-12)
    @test isapprox(accumulating.intercepted_snowfall, 0.003; rtol = 1e-12)
    @test abs(accumulating.snow_mass_flux) < 1e-12
    @test flooded.h > 0.1 && flooded.hs < 0.2
    @test isapprox(flooded.ice_mass_flux, -flooded.snow_mass_flux; rtol = 1e-12)
end
