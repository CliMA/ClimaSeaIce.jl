using ClimaSeaIce.SeaIceThermodynamics:
    BulkSalinityDiffusion,
    ColumnBasalHeatFlux,
    ColumnEnergyThermodynamics,
    SeaIceColumnDiscretization,
    ConductiveTemperatureTransport,
    FixedDrainedIceSalinityProfile,
    FixedSalinityBrinePocketEnergyRelation,
    FluxBoundary,
    FluxFunction,
    LinearizedSurfaceFlux,
    MeltingConstrainedFluxBalance,
    IceWaterThermalEquilibrium,
    MaykutUntersteinerConductivity,
    NoShortwaveAbsorption,
    PrescribedBulkSalinity,
    PrescribedTemperature,
    PrognosticBulkSalinity,
    QuadraticLiquidusEnergyRelation,
    column_energy_time_step!,
    column_height,
    complete_melt_energy,
    evolving_salinity_mushy_thermodynamics,
    internal_energy,
    melting_temperature,
    prescribed_salinity_enthalpy_thermodynamics,
    temperature
using ClimaSeaIce: ClimaSeaIce, SeaIceModel
using ClimaSeaIce: UniformThicknessDistribution, snow_slab_thermodynamics
using ClimaSeaIce.SeaIceThermodynamics: itd_factor
using Oceananigans
using Oceananigans.TimeSteppers: time_step!
using Oceananigans.Fields: interior, set!
using Oceananigans.Grids: znode
using Oceananigans: fields, prognostic_fields
using Test

@testset "Unconsolidated columns respect the basal reservoir and melting limit" begin
    for FT in (Float32, Float64), salinity in (0, 4),
        salinity_closure in (PrescribedBulkSalinity(FT(salinity)), PrognosticBulkSalinity()),
        initial_temperature in (-30, -10, -2), basal_reference in (-1.8, 0)
        grid = RectilinearGrid(FT; size = 4, z = SeaIceColumnDiscretization((0, 1)),
                               topology = (Flat, Flat, Bounded))
        relation = FixedSalinityBrinePocketEnergyRelation(FT)
        basal_temperature = FT(basal_reference)
        frazil_cooling = ifelse(basal_reference < 0, FT(20), FT(0))
        ice = ColumnEnergyThermodynamics(grid; relation, salinity_closure,
            salinity_transport = BulkSalinityDiffusion(FT; diffusivity = FT(1e-6)),
            energy_transport = ConductiveTemperatureTransport(FT; conductivity = FT(2)),
            heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = PrescribedTemperature(basal_temperature)))
        snow = snow_slab_thermodynamics(grid)
        model = SeaIceModel(grid; ice_thermodynamics = ice, snow_thermodynamics = snow,
            phase_transitions = relation.phase_transitions, sea_ice_density = FT(900),
            top_heat_flux = FT(0), bottom_heat_flux = ColumnBasalHeatFlux(FT(0), FT(0), -frazil_cooling),
            timestepper = :ForwardEuler)
        set!(model; h = 0, ℵ = 0, hs = 0)
        set!(ice; bulk_salinity = FT(salinity), temperature = FT(initial_temperature))
        Δt = FT(60)
        time_step!(model, Δt)
        volume = only(Array(interior(model.ice_thickness))) * only(Array(interior(model.ice_concentration)))
        latent_heat = ClimaSeaIce.SeaIceThermodynamics.latent_heat(relation.phase_transitions, basal_temperature)
        target_temperature = min(basal_temperature, melting_temperature(relation.phase_transitions.liquidus, FT(salinity)))
        @test volume ≈ frazil_cooling * Δt / (FT(900) * latent_heat)
        @test all(isapprox.(Array(interior(ice.fields.temperature)), target_temperature))
        @test all(isapprox.(Array(interior(ice.fields.bulk_salinity)), FT(salinity)))
        @test only(Array(interior(snow.top_surface_temperature))) ≈ target_temperature
    end
end

@testset "Float32 snow and ice surface balance at a 90-minute step" begin
    function thin_snow_step(FT)
        grid = RectilinearGrid(FT; size = (1, 1, 8), x = (0, 1), y = (0, 1),
                               z = SeaIceColumnDiscretization((0, 1)),
                               topology = (Bounded, Bounded, Bounded))
        relation = FixedSalinityBrinePocketEnergyRelation(FT)
        ice = prescribed_salinity_enthalpy_thermodynamics(grid;
            relation, salinity_profile = FixedDrainedIceSalinityProfile(FT),
            energy_transport = ConductiveTemperatureTransport(conductivity = MaykutUntersteinerConductivity(FT)),
            heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(),
                                        bottom = IceWaterThermalEquilibrium(salinity = 35.0)))
        snow = snow_slab_thermodynamics(grid)
        model = SeaIceModel(grid; ice_thermodynamics = ice, snow_thermodynamics = snow,
                            phase_transitions = relation.phase_transitions,
                            ice_consolidation_thickness = 0.05,
                            top_heat_flux = LinearizedSurfaceFlux(0.0, -1.0, 25.0),
                            bottom_heat_flux = 0.0, timestepper = :ForwardEuler)
        set!(model; h = 0.05, ℵ = 1.0, hs = 0.01, u = 0.0, v = 0.0)
        set!(ice; temperature = -1.0)
        time_step!(model, 5400.0)
        return (ice_temperature = first(interior(ice.auxiliary.surface_temperature)),
                snow_temperature = first(interior(snow.top_surface_temperature)),
                energy = Array(interior(ice.fields.internal_energy)))
    end

    f32 = thin_snow_step(Float32)
    f64 = thin_snow_step(Float64)
    @test all(isfinite, f32.energy)
    @test all(isfinite, (f32.ice_temperature, f32.snow_temperature))
    @test isapprox(f32.ice_temperature, f64.ice_temperature; atol = 0.02)
    @test isapprox(f32.snow_temperature, f64.snow_temperature; atol = 0.02)
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

    # Resolve the same ice/snow/atmosphere balance with one full 90-minute
    # thermodynamic step, without advancing the column in shorter pieces.
    time_step!(model, 90 * 60.0)
    Ts = snow.top_surface_temperature[1, 1, 1]
    Ti = ice.auxiliary.surface_temperature[1, 1, 1]
    hs = model.snow_thickness[1, 1, 1]
    @test all(isfinite, (Ts, Ti, hs))
    @test -30 < Ts < -20
    @test isapprox(0.31 * (Ti - Ts) / hs, atmosphere(1, 1, grid, Ts, model.clock, fields(model)); atol = 1e-5)
end

@testset "Column energy thermodynamics in SeaIceModel" begin
    grid = RectilinearGrid(size = 4,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary())

    evolving = evolving_salinity_mushy_thermodynamics(grid;
        relation,
        energy_transport = ConductiveTemperatureTransport(conductivity = 2.0, diffusivity = 1e-9),
        salinity_transport = BulkSalinityDiffusion(diffusivity = 2e-9),
        shortwave_absorption = NoShortwaveAbsorption(),
        heat_boundary_conditions)

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
    Tₘ = -relation.phase_transitions.liquidus.salinity_slope * S
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

@testset "Basal heat and frazil use grid-cell area" begin
    # One quarter of the cell has ice. The interface supplies 40 W/m² of ice,
    # hence only 10 W/m² of grid cell. Frazil is already per grid-cell area.
    grid = RectilinearGrid(size = 4, z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    ice = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation, salinity_profile = 4.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        heat_boundary_conditions = (top = FluxBoundary(), bottom = PrescribedTemperature(-5.0)))
    basal = ColumnBasalHeatFlux(40.0, 10.0, -20.0)
    model = SeaIceModel(grid; ice_thermodynamics = ice,
                        phase_transitions = relation.phase_transitions,
                        top_heat_flux = 0.0, bottom_heat_flux = basal)
    set!(model; h = 1.0, ℵ = 0.25)
    set!(ice; temperature = -5.0)

    Δt = 60.0
    volume_before = model.ice_thickness[1, 1, 1] * model.ice_concentration[1, 1, 1]
    time_step!(model, Δt)
    volume_after = model.ice_thickness[1, 1, 1] * model.ice_concentration[1, 1, 1]
    ρL = model.sea_ice_density[1, 1, 1] * relation.phase_transitions.reference_latent_heat

    # The +10 W/m² interface loss and -20 W/m² frazil gain yield net growth.
    @test isapprox(volume_after - volume_before, 10.0 * Δt / ρL; rtol = 1e-8)
end

@testset "set! syncs the column grid metric" begin
    grid = RectilinearGrid(size = 8, z = SeaIceColumnDiscretization((0, 2)), topology = (Flat, Flat, Bounded))
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 2.0),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = IceWaterThermalEquilibrium(salinity = 0)))
    model = SeaIceModel(grid; ice_thermodynamics = thermodynamics)

    set!(model, h = 0.3, ℵ = 1)
    @test isapprox(column_height(1, 1, grid), 0.3; rtol = 1e-12)
    copied = similar(grid)
    expanded = Oceananigans.Grids.with_halo((0, 0, 4), grid)
    expected_height = znode(1, 1, 4, grid, Center(), Center(), Center())
    @test znode(1, 1, 4, copied, Center(), Center(), Center()) ≈ expected_height
    @test znode(1, 1, 4, expanded, Center(), Center(), Center()) ≈ expected_height
end

@testset "Captured thin-snow surface temperature failure" begin
    # Arctic atmospheric state at model hour 1840.167.
    air_temperature = 277.728916422526
    air_humidity = 0.0048475598062699046
    pressure = 103132.7359375
    wind_speed = hypot(1.986268457770367, -0.03793264776467879)
    downward_longwave = 339.66806640625
    downward_shortwave = 57.385873349507946
    hs, h = 0.00023206826032621405, 1.0537630538084337
    Ttop, Tprevious = -0.2160028923486799, -0.216
    conductance = 2 * 0.1 / (h / 8)

    function atmospheric_flux(T)
        eice = 611.2 * exp(22.46 * T / (272.62 + T))
        qsat = 0.622 * eice / (pressure - 0.378 * eice)
        sensible = 1.3 * 1004.0 * 1.3e-3 * wind_speed * (T + 273.15 - air_temperature)
        latent = 1.3 * 2.83e6 * 1.3e-3 * wind_speed * (qsat - air_humidity)
        longwave = 0.97 * 5.670374419e-8 * (T + 273.15)^4 - downward_longwave
        shortwave = -0.2 * downward_shortwave
        return sensible + latent + longwave + shortwave
    end

    grid = RectilinearGrid(Float64; size = (1, 1, 8), x = (0, 1), y = (0, 1),
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Bounded, Bounded, Bounded))
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    ice = prescribed_salinity_enthalpy_thermodynamics(grid; relation, salinity_profile = 4,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.1),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = FluxBoundary()))
    snow = snow_slab_thermodynamics(grid)
    atmosphere = FluxFunction((i, j, grid, T, clock, fields) -> atmospheric_flux(T);
                              top_temperature_dependent = true)
    model = SeaIceModel(grid; ice_thermodynamics = ice, snow_thermodynamics = snow,
                        top_heat_flux = atmosphere, bottom_heat_flux = 0)
    set!(model; h, hs, ℵ = 1)
    set!(ice; temperature = Ttop)
    surface = (thickness = model.snow_thickness, density = model.snow_density,
               conductivity = 0.31, latent_heat = relation.phase_transitions.reference_latent_heat,
               itd_shape = nothing, temperature = snow.top_surface_temperature)
    column_energy_time_step!(ice, model.external_heat_fluxes, model.clock, fields(model),
                             model.ice_consolidation_thickness, 60; snow = surface)
    Ti = first(interior(ice.auxiliary.surface_temperature))
    Ts = first(interior(snow.top_surface_temperature))
    T = last(interior(ice.fields.temperature))
    Qmelt = first(interior(ice.auxiliary.surface_fluxes.snow_melt))
    Qsnow = atmospheric_flux(Ts) + Qmelt
    Qablation = first(interior(ice.auxiliary.surface_fluxes.ice_melt))
    @test -1 < Ti < 0
    @test abs(Qsnow + Qablation - conductance * (T - Ti)) < 1e-5
    @test abs(0.31 * (Ti - Ts) / hs - Qsnow) < 1e-5
    @test Qmelt == 0 && Qablation > 0

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
                       ice_conductivity = 2.0, ℵ⁰ = 1.0, bottom_heat_flux = 0.0, snow_conductivity = 0.31)
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
        snow = snow_slab_thermodynamics(grid; conductivity = snow_conductivity, itd_shape = snow_itd_shape)
        model = SeaIceModel(grid; ice_thermodynamics = ice, snow_thermodynamics = snow,
                            phase_transitions = relation.phase_transitions,
                            top_heat_flux = Qatm, bottom_heat_flux,
                            timestepper = :ForwardEuler, snowfall)
        set!(model; h = h⁰, ℵ = ℵ⁰, hs = hs⁰, u = 0.0, v = 0.0)
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
        skin_residual = hs > 1e-8 ? abs(snow_conductivity * itd_factor(snow_itd_shape, h⁰) * (Ti - Ts) / hs⁰ - Qatm) : 0.0
        skin_residual = Ts >= 0 ? 0.0 : skin_residual
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
    insulating = snow_case(-100.0, 0.2; snow_conductivity = 0)
    exhausted = snow_case(-2000.0, 1e-4; Δt = 20.0, initial_temperature = -0.3)
    fractional_growth = snow_case(-100.0, 0.01; Δt = 5400.0, initial_temperature = -0.3,
                                  h⁰ = 0.05, ℵ⁰ = 0.1,
                                  bottom_heat_flux = ColumnBasalHeatFlux(0.0, 0.0, -100.0))
    bare_melt = snow_case(-100.0, 1e-12; Δt = 900.0, h⁰ = 1.5, ice_conductivity = 0.0)
    accumulating = snow_case(0.0, 0.1; snowfall = 0.003)
    flooded = snow_case(0.0, 0.2; h⁰ = 0.1)
    cap_switch = snow_case(-500.0, 0.01; Δt = 5400.0, initial_temperature = -10.0, h⁰ = 0.15)
    snow_feedback = snow_case(-200.0, 0.01; Δt = 5400.0, initial_temperature = -10.0, h⁰ = 0.3)
    @test cold.hs == 0.2 && cold.Ts < 0
    @test distributed.Ts > cold.Ts
    @test isapprox(0.31 * UniformThicknessDistribution(Float64; categories = 4).conductivity_factor *
                   (distributed.Ti - distributed.Ts) / distributed.hs, 20.0; atol = 1e-5)
    @test distributed.skin_residual < 1e-6
    @test thin.hs ≈ 1e-7 && thin.Ts < 0 && thin.skin_residual < 1e-5
    @test dust.hs ≈ 1e-12 && isapprox(dust.Ts, dust.Ti; atol = 1e-8)
    @test warm.hs < 0.2 && warm.Ts == 0
    @test insulating.Ti ≈ -5
    @test isapprox(insulating.hs, 0.2 - 100 * 600 / (330 * 334000); atol = 1e-12)
    @test exhausted.hs < 1e-12 && isfinite(exhausted.Ts) && exhausted.ice_mass_flux < 0
    @test 0 < exhausted.ℵ < 1
    @test isapprox(exhausted.snow_mass_flux, -330 * 1e-4 / 20; rtol = 1e-12)
    phase_transitions = FixedSalinityBrinePocketEnergyRelation(Float64).phase_transitions
    Qmelt = 0.31 / 0.01 * fractional_growth.Ti + 100
    @test fractional_growth.ℵ > 0.1
    @test fractional_growth.Ts == 0
    @test 0 < Qmelt < 330 * phase_transitions.reference_latent_heat * 0.01 / 5400
    @test isapprox(fractional_growth.snow_mass_flux,
                   -0.1 * Qmelt / phase_transitions.reference_latent_heat; rtol = 1e-10)
    Tₘ = -phase_transitions.liquidus.salinity_slope * 4.0
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
    @test all(isfinite, (cap_switch.h, cap_switch.hs, cap_switch.Ti,
                         snow_feedback.h, snow_feedback.hs, snow_feedback.Ti))
    @test cap_switch.h > 0 && cap_switch.hs < 0.01
    @test snow_feedback.h > 0 && snow_feedback.hs < 0.01 && snow_feedback.Ti < 0
    @test max(cap_switch.skin_residual, snow_feedback.skin_residual) < 1e-6
end

@testset "One-cell snow melt has an analytical heat and mass budget" begin
    grid = RectilinearGrid(size = 1, z = SeaIceColumnDiscretization((0, 1)), topology = (Flat, Flat, Bounded))
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    ice = prescribed_salinity_enthalpy_thermodynamics(grid; relation, salinity_profile = 0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 5),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = FluxBoundary()))
    snow = snow_slab_thermodynamics(grid)
    model = SeaIceModel(grid; ice_thermodynamics = ice, snow_thermodynamics = snow,
                        top_heat_flux = -100, bottom_heat_flux = 0, phase_transitions = relation.phase_transitions)
    set!(model; h = 1, ℵ = 1, hs = 0.01)
    set!(ice; temperature = -8.2)
    Δt = 600
    ρc = relation.phase_transitions.density * relation.phase_transitions.heat_capacity
    C = 10
    R = 0.01 / 0.31
    B = C / (1 + C * R)
    T = -8.2 / (1 + Δt * B / ρc)
    Ti = T * C * R / (1 + C * R)
    Qs = 100 + 0.31 * Ti / 0.01
    hs = 0.01 - Δt * Qs / (330 * relation.phase_transitions.reference_latent_heat)
    time_step!(model, Δt)
    @test first(interior(snow.top_surface_temperature)) == 0
    @test isapprox(first(interior(ice.fields.temperature)), T; atol = 1e-10)
    @test isapprox(first(interior(ice.auxiliary.surface_temperature)), Ti; atol = 1e-10)
    @test isapprox(first(interior(model.snow_thickness)), hs; atol = 1e-12)
    @test isapprox(ρc * (T + 8.2) + 330 * relation.phase_transitions.reference_latent_heat * (0.01 - hs), 100Δt; atol = 1e-7)

    set!(model; h = 0.04, ℵ = 0.5, hs = 0.01)
    set!(ice; temperature = 0)
    time_step!(model, Δt)
    @test isapprox(first(interior(model.ice_thickness)) * first(interior(model.ice_concentration)), 0.02; atol = 1e-12)
    @test isapprox(first(interior(model.snow_thickness)), 0.01 - 100Δt / (330 * relation.phase_transitions.reference_latent_heat); atol = 1e-12)
    @test abs(first(interior(model.mass_fluxes.thermodynamics.ice))) < 1e-12
end

@testset "A thin snow-covered column from a failed global step converges" begin
    # Column (292, 47) at iteration 905 of a global run whose solver backtracked on the normalized defect, and failed here:
    # along its second step the residual falls at every fraction, but the defect rises.
    E⁰ = [25723765.02132106, 27658415.953812763, 29652550.64854156, 31579968.399667785,
          32999238.126846775, 32730876.60558527, 28118057.369986087, 14615218.968743838]
    h, hs, Sₒ, Δt = 0.050359075912101384, 0.0035281146699675274, 34.157656051293074, 5400

    grid = RectilinearGrid(Float64; size = (1, 1, 8), x = (0, 1), y = (0, 1),
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Bounded, Bounded, Bounded))
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    ice = prescribed_salinity_enthalpy_thermodynamics(grid; relation,
        salinity_profile = FixedDrainedIceSalinityProfile(Float64),
        energy_transport = ConductiveTemperatureTransport(conductivity = MaykutUntersteinerConductivity(Float64)),
        heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = IceWaterThermalEquilibrium(Sₒ)))
    snow = snow_slab_thermodynamics(grid)
    atmosphere = LinearizedSurfaceFlux(-60.93563649479057, -0.01198902224905396, 40.5078141969176)
    model = SeaIceModel(grid; ice_thermodynamics = ice, snow_thermodynamics = snow,
                        top_heat_flux = atmosphere, bottom_heat_flux = 0)
    set!(model; h, hs, ℵ = 1)
    set!(ice; internal_energy = reshape(E⁰, 1, 1, 8))
    surface = (thickness = model.snow_thickness, density = model.snow_density,
               conductivity = 0.31, latent_heat = relation.phase_transitions.reference_latent_heat,
               itd_shape = nothing, temperature = snow.top_surface_temperature)

    column_energy_time_step!(ice, model.external_heat_fluxes, model.clock, fields(model),
                             model.ice_consolidation_thickness, Δt; snow = surface)

    E = vec(interior(ice.fields.internal_energy))
    T = vec(interior(ice.fields.temperature))
    zᵇ = -(900h + 330hs) / relation.phase_transitions.liquid_density
    Tᵇ = melting_temperature(relation.phase_transitions.liquidus, Sₒ, zᵇ)
    Δz = h / 8
    Qᵇ = 2 * first(interior(ice.auxiliary.thermal_conductivity)) / Δz * (first(T) - Tᵇ)
    Qᶜ = first(interior(ice.auxiliary.surface_fluxes.heat))

    @test abs(Δz * sum(E .- E⁰) / Δt + Qᶜ + Qᵇ) < 1e-3
    @test all(T .< 0)
end

@testset "Unconsolidated snow and ice share the basal-temperature reservoir" begin
    for initial_temperature in (-30.0, -10.0, -2.0)
        grid = RectilinearGrid(size = 8, z = SeaIceColumnDiscretization((0, 1)), topology = (Flat, Flat, Bounded))
        relation = FixedSalinityBrinePocketEnergyRelation(Float64)
        Tb = -relation.phase_transitions.liquidus.salinity_slope * 34
        ice = prescribed_salinity_enthalpy_thermodynamics(grid; relation, salinity_profile = 0,
            energy_transport = ConductiveTemperatureTransport(conductivity = 2),
            heat_boundary_conditions = (top = MeltingConstrainedFluxBalance(), bottom = PrescribedTemperature(Tb)))
        model = SeaIceModel(grid; ice_thermodynamics = ice, snow_thermodynamics = snow_slab_thermodynamics(grid),
                            phase_transitions = relation.phase_transitions, sea_ice_density = 917,
                            top_heat_flux = -200, bottom_heat_flux = 0, timestepper = :ForwardEuler)
        set!(model; h = 0.04, ℵ = 1, hs = 0.005)
        set!(ice; temperature = initial_temperature)
        Δt = 60
        Qsnow = 200 + 0.31 * Tb / 0.005
        time_step!(model, Δt)
        V = first(interior(model.ice_thickness)) * first(interior(model.ice_concentration))
        Vs = first(interior(model.snow_thickness)) * first(interior(model.ice_concentration))
        snow_heat = 330 * relation.phase_transitions.reference_latent_heat * (0.005 - Vs) / Δt
        @test isapprox(first(interior(ice.auxiliary.surface_temperature)), Tb; atol = 1e-12)
        @test isapprox(snow_heat, Qsnow; atol = 1e-6)
        @test isapprox(917 * ClimaSeaIce.SeaIceThermodynamics.latent_heat(relation.phase_transitions, Tb) * (0.04 - V) / Δt + snow_heat, 200; atol = 1e-6)
    end
end

@testset "Column ice with slab snow is validated when the model is built" begin
    grid = RectilinearGrid(size = 8, z = SeaIceColumnDiscretization((0, 1)), topology = (Flat, Flat, Bounded))
    ice(top) = prescribed_salinity_enthalpy_thermodynamics(grid; salinity_profile = 0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 2),
        heat_boundary_conditions = (top, bottom = PrescribedTemperature(-1.8)))
    snow = snow_slab_thermodynamics(grid)
    @test SeaIceModel(grid; ice_thermodynamics = ice(MeltingConstrainedFluxBalance()), snow_thermodynamics = snow) isa SeaIceModel
    @test_throws ArgumentError SeaIceModel(grid; ice_thermodynamics = ice(PrescribedTemperature(-5)), snow_thermodynamics = snow)
end
