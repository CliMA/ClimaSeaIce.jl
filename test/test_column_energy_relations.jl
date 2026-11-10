using ClimaSeaIce.SeaIceThermodynamics:
    BubblyBrineConductivity,
    FixedDrainedIceSalinityProfile,
    FixedSalinityBrinePocketEnergyRelation,
    MaykutUntersteinerConductivity,
    QuadraticLiquidusEnergyRelation,
    brine_salinity,
    complete_melt_energy,
    face_thermal_conductivity,
    ice_thermal_conductivity,
    internal_energy,
    liquid_fraction,
    salinity_at_normalized_depth,
    salinity_at_normalized_height,
    temperature,
    temperature_energy_derivative,
    temperature_salinity_derivative
using Oceananigans
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

        @test test_liquidus_consistency(Float32) < 2e-6
        test_liquid_fraction_bounds(Float32)
        @test test_complete_melt_energy(Float32) < 1e-6
    end
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
end
