using Oceananigans
using ClimaSeaIce

include("distributed_tests_utils.jl")

test_distributed_halos = """
    using MPI
    MPI.Init()

    include("distributed_tests_utils.jl")
    test_extended_halos()
"""

run_distributed_grid = """
    using MPI
    MPI.Init()
    include("distributed_tests_utils.jl")

    test_distributed_simulations()
"""

run_distributed_jld2writer = """
    using MPI
    MPI.Init()

    include("distributed_tests_utils.jl")
    arch = Distributed(CPU(), partition = Partition(2, 2))
    run_distributed_sea_ice_jl2dwriter(arch, "distributed_jld2_writer.jld2")
"""

@testset "Test distributed sea-ice grid simulations..." begin

    # @info "Testing (4, 1) Distributed sea ice halo extension"

    # # Test extended halos
    # write("distributed_halo_tests.jl", test_distributed_halos)
    # run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project -O0 distributed_halo_tests.jl`)
    # rm("distributed_halo_tests.jl")

    # # Run the serial computation
    # @info "Testing 4-ranks Distributed sea ice simulations equality"

    # # Run the distributed grid simulation 
    # write("distributed_simulation_tests.jl", run_distributed_grid)
    # run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project -O0 distributed_simulation_tests.jl`)
    # rm("distributed_simulation_tests.jl")

    # @info "Testing JLD2Writer on distributed sea ice simulations"

    write("distributed_jld2writer_tests.jl", run_distributed_jld2writer)
    run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project -O0 distributed_jld2writer_tests.jl`)
    rm("distributed_jld2writer_tests.jl")

    @test isfile("distributed_jld2_writer_rank0.jld2")
    @test isfile("distributed_jld2_writer_rank1.jld2")
    @test isfile("distributed_jld2_writer_rank2.jld2")
    @test isfile("distributed_jld2_writer_rank3.jld2")

    # This does not work, see https://github.com/CliMA/Oceananigans.jl/pull/5188 
    h = FieldTimeSeries("distributed_jld2_writer_rank0.jld2", "h")
    ℵ = FieldTimeSeries("distributed_jld2_writer_rank0.jld2", "ℵ")

    @test length(h.times) == 21
    @test length(ℵ.times) == 21

    rm("distributed_jld2_writer_rank0.jld2")
    rm("distributed_jld2_writer_rank1.jld2")
    rm("distributed_jld2_writer_rank2.jld2")
    rm("distributed_jld2_writer_rank3.jld2")
end
