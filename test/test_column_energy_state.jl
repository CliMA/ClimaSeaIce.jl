using ClimaSeaIce.SeaIceThermodynamics:
    BulkSalinityDiffusion,
    ColumnEnergyThermodynamics,
    SeaIceColumnDiscretization,
    ConductiveTemperatureTransport,
    ExponentialShortwaveAbsorption,
    FixedDrainedIceSalinityProfile,
    FluxBoundary,
    NoSalinityTransport,
    NoShortwaveAbsorption,
    PrescribedBulkSalinity,
    PrognosticBulkSalinity,
    QuadraticLiquidusEnergyRelation,
    brine_salinity,
    evolving_salinity_mushy_thermodynamics,
    initialize_column_interfaces!,
    internal_energy,
    liquid_fraction,
    prescribed_salinity_enthalpy_thermodynamics,
    temperature,
    temperature_energy_derivative,
    temperature_salinity_derivative
using Adapt
using Oceananigans
using Oceananigans.Fields: interior, set!
using Oceananigans: fields, prognostic_fields, prognostic_state, restore_prognostic_state!
using Test

@testset "Fixed drained-ice salinity profile on a column grid" begin
    profile = FixedDrainedIceSalinityProfile(Float64)

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
end
