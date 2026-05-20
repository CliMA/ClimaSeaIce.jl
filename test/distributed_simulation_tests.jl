    using MPI
    MPI.Init()
    include("distributed_tests_utils.jl")

    test_distributed_simulations()
