using ClimaSeaIce
using Oceananigans
using Oceananigans: prognostic_fields
using Oceananigans.Architectures: architecture

@testset "Sea ice advection" begin
    @info "Running sea ice advection test"

    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))

    model = SeaIceModel(grid, advection=WENO())

    @test architecture(model) == architecture(grid)

    @test !(model.velocities.u isa Nothing)
    @test !(model.velocities.v isa Nothing)

    # test that model runs
    @test begin
        time_step!(model, 1)
        true
    end
end

@testset "Sea ice advection with snow" begin
    @info "Running sea ice advection with snow test"

    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))
    snow_thermo = SlabSnowThermodynamics(grid)
    model = SeaIceModel(grid; advection=WENO(), snow_thermodynamics=snow_thermo)

    set!(model, h=1, ℵ=1, hs=0.2)

    @test model.snow_thickness isa Field
    @test :hs in keys(prognostic_fields(model))

    @test begin
        time_step!(model, 1)
        true
    end
end

@testset "Sea ice momentum equations" begin
    @info "Running sea ice momentum equations test"

    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))
    dynamics = SeaIceMomentumEquation(grid, rheology=ViscousRheology(ν=1000))

    model = SeaIceModel(grid; dynamics, ice_thermodynamics=nothing, advection=WENO())

    @test !(model.velocities.u isa Nothing)
    @test !(model.velocities.v isa Nothing)

    # test that model runs
    @test begin
        time_step!(model, 1)
        true
    end
end
