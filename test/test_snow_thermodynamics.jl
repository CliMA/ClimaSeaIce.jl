using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: ConductiveFlux, PhaseTransitions, GammaThicknessDistribution,
    IceSnowConductiveFlux, ice_snow_conductive_flux, interface_temperature, latent_heat, slab_internal_heat_flux, itd_factor,
    bottom_temperature, melting_temperature
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: PrescribedTemperature, FluxFunction, IceWaterThermalEquilibrium
using Oceananigans
using Oceananigans: prognostic_fields
using Oceananigans.Fields: interior
using Test

@testset "Snow model construction" begin
    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))

    # Model with snow
    snow_thermo = snow_slab_thermodynamics(grid)
    model = SeaIceModel(grid; snow_thermodynamics=snow_thermo)
    @test model.snow_thermodynamics isa SlabThermodynamics
    @test model.snow_thickness isa Field

    # Both layers store the raw conductive-flux coefficient; the combined
    # snow+ice flux (IceSnowConductiveFlux) is assembled inline in the
    # layered kernel, not stored on the thermodynamics.
    snow_flux = model.snow_thermodynamics.internal_heat_flux
    @test snow_flux isa ConductiveFlux
    @test snow_flux.conductivity ≈ 0.31

    ice_flux = model.ice_thermodynamics.internal_heat_flux
    @test ice_flux isa ConductiveFlux

    # Model without snow
    model_no_snow = SeaIceModel(grid)
    @test isnothing(model_no_snow.snow_thermodynamics)
    @test isnothing(model_no_snow.snow_thickness)
end

@testset "Backward compatibility without snow" begin
    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))
    model = SeaIceModel(grid)
    @test isnothing(model.snow_thermodynamics)
    @test isnothing(model.snow_thickness)

    set!(model, h=1, ℵ=1)
    simulation = Simulation(model, Δt=1.0, stop_iteration=3)
    run!(simulation)
    @test model.clock.iteration == 3
end

@testset "Snow insulation effect" begin
    ki = 2.0
    ks = 0.31
    hi = 1.0
    hs = 0.3
    Tu = -10.0
    Tb = -1.8

    # Ice-only conductive flux
    Fc_no_snow = -ki * (Tu - Tb) / hi

    # Combined snow+ice conductive flux (resistors in series)
    R = hs / ks + hi / ki
    Fc_with_snow = (Tb - Tu) / R

    # Snow should reduce the magnitude of the conductive flux
    @test abs(Fc_with_snow) < abs(Fc_no_snow)

    # With zero snow, should recover the no-snow result
    R_zero = 0.0 / ks + hi / ki
    Fc_zero_snow = (Tb - Tu) / R_zero
    @test Fc_zero_snow ≈ Fc_no_snow

    # Thicker snow -> even less flux
    R_thick = 1.0 / ks + hi / ki
    Fc_thick_snow = (Tb - Tu) / R_thick
    @test abs(Fc_thick_snow) < abs(Fc_with_snow)
end

@testset "Interface temperature" begin
    ki = 2.0
    ks = 0.31
    hi = 1.0
    hs = 0.3
    Tu = -10.0
    Tb = -1.8

    Ri = hi / ki
    Rs = hs / ks
    R  = Rs + Ri

    Tsi = Tb + (Tu - Tb) * Ri / R

    # Interface temperature should be between Tu and Tb
    @test Tsi > Tu
    @test Tsi < Tb

    # With no snow (hs = 0): Tsi = Tu
    Tsi_no_snow = Tb + (Tu - Tb) * Ri / Ri  # R = Ri when hs = 0
    @test Tsi_no_snow ≈ Tu
end

@testset "Snow-ice formation (flooding)" begin
    grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

    snow_thermo = snow_slab_thermodynamics(grid)
    ice_thermo  = SlabThermodynamics(grid; top_heat_boundary_condition = PrescribedTemperature(-5.0))

    model = SeaIceModel(grid;
                        ice_thermodynamics = ice_thermo,
                        snow_thermodynamics = snow_thermo)

    # Heavy snow on thin ice -> negative freeboard
    hi = 0.5
    hs = 1.0
    set!(model, h=hi, ℵ=1, hs=hs)

    time_step!(model, 1)

    hi⁺ = first(interior(model.ice_thickness))
    hs⁺ = first(interior(model.snow_thickness))

    # After flooding, ice increases and snow decreases
    @test hi⁺ > hi
    @test hs⁺ < hs
end

@testset "Snowfall accumulation" begin
    grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

    snow_thermo = snow_slab_thermodynamics(grid)
    ice_thermo  = SlabThermodynamics(grid)

    Ps = 1e-5  # kg/m²/s snowfall rate
    model = SeaIceModel(grid;
                        ice_thermodynamics = ice_thermo,
                        snow_thermodynamics = snow_thermo,
                        snowfall = Ps)

    set!(model, h=1, ℵ=1, hs=0)

    Δt = 3600  # 1 hour
    time_step!(model, Δt)

    hs⁺ = first(interior(model.snow_thickness))
    # Snow should have accumulated
    @test hs⁺ > 0
end

@testset "Snow melts before ice" begin
    grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

    snow_thermo = snow_slab_thermodynamics(grid)
    ice_thermo  = SlabThermodynamics(grid)

    # Negative top_heat_flux means incoming heat (solar radiation), which
    # drives the surface to the melting point and creates a flux imbalance.
    model = SeaIceModel(grid;
                        ice_thermodynamics = ice_thermo,
                        snow_thermodynamics = snow_thermo,
                        top_heat_flux = -100) # W/m² incoming

    hi = 2.0
    hs = 0.1
    set!(model, h=hi, ℵ=1, hs=hs)

    time_step!(model, 3600)  # 1 hour

    hs⁺ = first(interior(model.snow_thickness))

    # Snow should decrease (absorbs top melting energy first)
    @test hs⁺ < hs
end

@testset "Time stepping with snow" begin
    for timestepper in (:ForwardEuler, :SplitRungeKutta3)
        grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))

        snow_thermo = snow_slab_thermodynamics(grid)

        model = SeaIceModel(grid;
            snow_thermodynamics = snow_thermo,
            snow_density = 330,
            advection = WENO(),
            timestepper)

        set!(model, h=1, ℵ=1, hs=0.1)
        simulation = Simulation(model, Δt=1.0, stop_iteration=3)
        run!(simulation)
        @test model.clock.iteration == 3
    end
end

@testset "Sub-grid thickness correction of the conductivity" begin
    # N equal-area categories of thickness (2i - 1) h / N conduct Σᵢ 1 / (2i - 1) times more than the mean thickness
    bare = ConductiveFlux(Float64; conductivity=2)
    flux = ConductiveFlux(Float64; conductivity=2, itd_shape=UniformThicknessDistribution(Float64; categories=5))
    @test slab_internal_heat_flux(bare, -10.0, -1.8, 1.0) ≈ 2 * 8.2
    @test slab_internal_heat_flux(flux, -10.0, -1.8, 1.0) ≈ 2 * 8.2 * (1 + 1/3 + 1/5 + 1/7 + 1/9)
    @test UniformThicknessDistribution(Float64; categories=1).conductivity_factor == 1

    # The thickness-dependent correction s / (s - 1) grows as the assumed distribution widens with thickness
    flux = ConductiveFlux(Float64; conductivity=2, itd_shape=GammaThicknessDistribution())
    thin_conductance  = 0.1 * slab_internal_heat_flux(flux, -10.0, -1.8, 0.1)
    thick_conductance = 5.0 * slab_internal_heat_flux(flux, -10.0, -1.8, 5.0)
    @test thick_conductance > thin_conductance > 2 * 8.2
end

@testset "Sub-grid thickness correction of the snow and ice column" begin
    grid = RectilinearGrid(size=(1, 1), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))
    model = SeaIceModel(grid; snow_thermodynamics=snow_slab_thermodynamics(grid))
    set!(model, h=1, hs=0.3)

    fields = (h = model.ice_thickness, hs = model.snow_thickness, ρi = model.sea_ice_density, ρs = model.snow_density)
    liquidus = model.phase_transitions.liquidus
    bottom_heat_boundary_condition = model.ice_thermodynamics.heat_boundary_conditions.bottom
    column_flux(flux) = ice_snow_conductive_flux(1, 1, grid, -10.0, model.clock, fields, (; flux, liquidus, bottom_heat_boundary_condition))
    snow_ice_temperature(flux) = interface_temperature(1, 1, grid, flux, bottom_heat_boundary_condition, liquidus, -10.0, fields)

    # Snow and ice are scaled alike: the column conducts `itd_factor` times more and the snow-ice interface stays put
    bare = IceSnowConductiveFlux(0.31, 2.0)
    for itd_shape in (UniformThicknessDistribution(Float64; categories=5), GammaThicknessDistribution(Float64))
        corrected = IceSnowConductiveFlux(0.31, 2.0, itd_shape)
        @test column_flux(corrected) ≈ itd_factor(itd_shape, 1.0) * column_flux(bare)
        @test snow_ice_temperature(corrected) ≈ snow_ice_temperature(bare)
    end
end

@testset "Melting temperature at the base of floating ice" begin
    grid = RectilinearGrid(size=(1, 1), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))
    bottom_heat_boundary_condition = IceWaterThermalEquilibrium(Float64; salinity=34)
    ice_thermodynamics = sea_ice_slab_thermodynamics(grid; bottom_heat_boundary_condition)

    bare_ice = SeaIceModel(grid; ice_thermodynamics, sea_ice_density=900)
    snowy_ice = SeaIceModel(grid; ice_thermodynamics, sea_ice_density=900, snow_density=330,
                            snow_thermodynamics=snow_slab_thermodynamics(grid))

    set!(bare_ice, h=2)
    set!(snowy_ice, h=2, hs=0.3)

    liquidus = bare_ice.phase_transitions.liquidus
    Tb(model) = bottom_temperature(1, 1, grid, bottom_heat_boundary_condition, liquidus, Oceananigans.fields(model))

    # The base of floating ice sits at the depth where the water carries the weight of the ice and snow above it
    @test Tb(bare_ice)  ≈ melting_temperature(liquidus, 34, - 900 * 2 / 1026)
    @test Tb(snowy_ice) ≈ melting_temperature(liquidus, 34, - (900 * 2 + 330 * 0.3) / 1026)
    @test Tb(snowy_ice) < Tb(bare_ice) < melting_temperature(liquidus, 34)

    set!(bare_ice, h=0)
    @test Tb(bare_ice) == melting_temperature(liquidus, 34)
end
