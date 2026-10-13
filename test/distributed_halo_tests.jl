    using MPI
    MPI.Init()
    include("distributed_tests_utils.jl")
    test_extended_halos()
