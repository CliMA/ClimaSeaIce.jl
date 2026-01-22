using Oceananigans
using ClimaSeaIce

include("distributed_tests_utils.jl")

run_slab_distributed_grid = """
    using MPI
    MPI.Init()

    include("distributed_tests_utils.jl")
    arch = Distributed(CPU(), partition = Partition(1, 4))
    run_distributed_sea_ice(arch, "distributed_yslab_seaice.jld2")
"""

run_pencil_distributed_grid = """
    using MPI
    MPI.Init()

    include("distributed_tests_utils.jl")
    arch = Distributed(CPU(), partition = Partition(2, 2))
    run_distributed_sea_ice(arch, "distributed_pencil_seaice.jld2")
"""

run_large_pencil_distributed_grid = """
    using MPI
    MPI.Init()

    include("distributed_tests_utils.jl")
    arch = Distributed(CPU(), partition = Partition(4, 2))
    run_distributed_sea_ice(arch, "distributed_large_pencil_seaice.jld2")
"""

run_distributed_jld2writer = """
    using MPI
    MPI.Init()

    include("distributed_tests_utils.jl")
    arch = Distributed(CPU(), partition = Partition(2, 2))
    run_distributed_sea_ice_jl2dwriter(arch, "distributed_jld2_writer.jld2")
"""

@testset "Test distributed sea-ice grid simulations..." begin
    # Run the serial computation
    grid = RectilinearGrid(CPU(); 
                           size = (100, 100, 1), 
                           x = (-10kilometers, 10kilometers), 
                           y = (-10kilometers, 10kilometers), 
                           z = (-1, 0), 
                           halo = (5, 5, 5))

    model = run_distributed_simulation(grid)

    # Retrieve Serial quantities
    us, vs = model.velocities
    hs = model.ice_thickness
    ℵs = model.ice_concentration

    us = interior(us, :, :, 1)
    vs = interior(vs, :, :, 1)
    hs = interior(hs, :, :, 1)
    ℵs = interior(ℵs, :, :, 1)

    @info "Testing (4, 1) Distributed sea ice simulations equality"

    # Run the distributed grid simulation with a slab configuration
    write("distributed_slab_tests.jl", run_slab_distributed_grid)
    run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project -O0 distributed_slab_tests.jl`)
    rm("distributed_slab_tests.jl")

    # Retrieve Parallel quantities
    up = jldopen("distributed_yslab_seaice.jld2")["u"]
    vp = jldopen("distributed_yslab_seaice.jld2")["v"]
    hp = jldopen("distributed_yslab_seaice.jld2")["h"]
    ℵp = jldopen("distributed_yslab_seaice.jld2")["ℵ"]

    rm("distributed_yslab_seaice.jld2")

    # Test slab partitioning
    @test all(us .≈ up)
    @test all(vs .≈ vp)
    @test all(hs .≈ hp)
    @test all(ℵs .≈ ℵp)

    @info "Testing (2, 2) Distributed sea ice simulations equality"

    # Run the distributed grid simulation with a pencil configuration
    write("distributed_tests.jl", run_pencil_distributed_grid)
    run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project -O0 distributed_tests.jl`)
    rm("distributed_tests.jl")

    # Retrieve Parallel quantities
    up = jldopen("distributed_pencil_seaice.jld2")["u"]
    vp = jldopen("distributed_pencil_seaice.jld2")["v"]
    hp = jldopen("distributed_pencil_seaice.jld2")["h"]
    ℵp = jldopen("distributed_pencil_seaice.jld2")["ℵ"]

    rm("distributed_pencil_seaice.jld2")

    @test all(us .≈ up)
    @test all(vs .≈ vp)
    @test all(hs .≈ hp)
    @test all(ℵs .≈ ℵp)

    @info "Testing (4, 2) Distributed sea ice simulations equality"

    write("distributed_large_pencil_tests.jl", run_large_pencil_distributed_grid)
    run(`$(mpiexec()) -n 8 $(Base.julia_cmd()) --project -O0 distributed_large_pencil_tests.jl`)
    rm("distributed_large_pencil_tests.jl")

    # Retrieve Parallel quantities
    up = jldopen("distributed_large_pencil_seaice.jld2")["u"]
    vp = jldopen("distributed_large_pencil_seaice.jld2")["v"]
    hp = jldopen("distributed_large_pencil_seaice.jld2")["h"]
    ℵp = jldopen("distributed_large_pencil_seaice.jld2")["ℵ"]

    rm("distributed_large_pencil_seaice.jld2")

    @test all(us .≈ up)
    @test all(vs .≈ vp)
    @test all(hs .≈ hp)
    @test all(ℵs .≈ ℵp)

    @info "Testing JLD2Writer on distributed sea ice simulations"

    write("distributed_jld2writer_tests.jl", run_distributed_jld2writer)
    run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project -O0 distributed_jld2writer_tests.jl`)
    rm("distributed_jld2writer_tests.jl")

    @test isfile("distributed_jld2_writer_rank0.jld2")
    @test isfile("distributed_jld2_writer_rank1.jld2")
    @test isfile("distributed_jld2_writer_rank2.jld2")
    @test isfile("distributed_jld2_writer_rank3.jld2")

    # This does not work, see https://github.com/CliMA/Oceananigans.jl/pull/5188 
    # h = FieldTimeSeries("distributed_jld2_writer.jld2", "h")
    # ℵ = FieldTimeSeries("distributed_jld2_writer.jld2", "ℵ")
    h = FieldTimeSeries("distributed_jld2_writer_rank0.jld2", "h")
    ℵ = FieldTimeSeries("distributed_jld2_writer_rank0.jld2", "ℵ")

    @test length(h.times) == 20
    @test length(ℵ.times) == 20

    rm("distributed_jld2_writer_rank0.jld2")
    rm("distributed_jld2_writer_rank1.jld2")
    rm("distributed_jld2_writer_rank2.jld2")
    rm("distributed_jld2_writer_rank3.jld2")
end
