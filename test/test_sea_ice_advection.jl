using ClimaSeaIce
using Oceananigans

@testset "Sea ice advection" begin
    @info "Running sea ice advection test"

    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))

    model = SeaIceModel(grid, advection=WENO()) 

    @test !(model.velocities.u isa Nothing)
    @test !(model.velocities.v isa Nothing)

    # test that model runs with RK3
    @test begin
        time_step!(model, 1)
        true
    end

    model = SeaIceModel(grid, advection=WENO(), timestepper=:QuasiAdamsBashforth2) 

    @test !(model.velocities.u isa Nothing)
    @test !(model.velocities.v isa Nothing)

    # test that model runs with AB2
    @test begin
        time_step!(model, 1)
        true
    end
end
