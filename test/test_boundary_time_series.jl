using Test
using Oceananigans
using ClimaSeaIce

using Oceananigans.Architectures: on_architecture
using Oceananigans.OutputReaders: InMemory, OnDisk

@testset "Boundary time series advance with the sea ice clock" begin
    grid = RectilinearGrid(size=(4, 4), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))

    times = collect(0.0:1.0:5.0)
    path = joinpath(mktempdir(), "boundary_thickness.jld2")
    ondisk = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; backend=OnDisk(), path, name="h")
    sample = Field{Nothing, Nothing, Nothing}(grid)

    for n in eachindex(times)
        set!(sample, n)
        set!(ondisk, sample, n)
    end

    # Only two frames are resident, so the window must move as the clock does
    h_boundary = FieldTimeSeries(path, "h"; backend=InMemory(2))
    h_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(h_boundary))

    model = SeaIceModel(grid; dynamics=nothing, ice_thermodynamics=nothing, advection=nothing,
                        boundary_conditions=(; h=h_bcs))

    Nx = size(grid, 1)

    for n in 1:4
        time_step!(model, 1)
        h = on_architecture(CPU(), model.ice_thickness)
        # The interior stays at zero, so the halo is twice the boundary value, which is time + 1
        @test h[Nx+1, 2, 1] ≈ 2 * (model.clock.time + 1)
    end
end
