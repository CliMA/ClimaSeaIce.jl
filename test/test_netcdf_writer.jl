using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics
using ClimaSeaIce.SeaIceThermodynamics
using Test

using NCDatasets
using Oceananigans
using Oceananigans.OutputWriters: NetCDFWriter

function test_netcdf_writer(arch)
    Nx, Ny = 16, 16
    Lx, Ly = 100, 100
    Δt = 1

    grid = RectilinearGrid(arch, size=(Nx, Ny), x=(0, Lx), y=(0, Ly), topology=(Bounded, Bounded, Flat))

    for dynamics in (nothing, SeaIceMomentumEquation(grid))
        ice_thermodynamics = SlabSeaIceThermodynamics(grid)
        model = SeaIceModel(grid; dynamics, ice_thermodynamics)

        set!(model, h=0.5, ℵ=1)

        simulation = Simulation(model, Δt=Δt, stop_iteration=5)

        filepath = "test_sea_ice_netcdf.nc"

        simulation.output_writers[:nc] = NetCDFWriter(model, Oceananigans.prognostic_fields(model);
                                                      filename = filepath,
                                                      schedule = IterationInterval(1),
                                                      overwrite_existing = true)

        run!(simulation)

        ds = NCDataset(filepath)

        # Check that the time dimension exists
        @test haskey(ds.dim, "time")
        @test length(ds["time"]) == 6  # iteration 0 through 5

        # Check that sea ice fields are present
        @test haskey(ds, "h")
        @test haskey(ds, "ℵ")

        # Check that velocity fields are present when dynamics are enabled
        if !isnothing(dynamics)
            @test haskey(ds, "u")
            @test haskey(ds, "v")
        end

        close(ds)
        rm(filepath, force=true)
    end
end

@testset "NetCDFWriter Tests" begin
    test_netcdf_writer(CPU())
end
