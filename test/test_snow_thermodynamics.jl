using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: ConductiveFlux, PhaseTransitions,
    IceSnowConductiveFlux, ice_snow_conductive_flux, interface_temperature, latent_heat
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: PrescribedTemperature, FluxFunction
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
