using MPI
using JLD2 

distributed_sea_ice = """
    using ClimaSeaIce
    using ClimaSeaIce.SeaIceThermodynamics
    using ClimaSeaIce.SeaIceDynamics
    using ClimaSeaIce.Rheologies
    using Oceananigans
    using Oceananigans.DistributedComputations
    using JLD2

    arch = Distributed(CPU(); partition=Partition(x=2, y=2))

    grid_2d = RectilinearGrid(arch, size=(10, 10),     x=(0, 1), y=(0, 1),            topology=(Bounded, Bounded, Flat))
    grid_3d = RectilinearGrid(arch, size=(10, 10, 10), x=(0, 1), y=(0, 1), z=(-1, 0), topology=(Bounded, Bounded, Bounded))

    iterations = Int[]

    for grid in (grid_2d, grid_3d)
        rheologies = (ElastoViscoPlasticRheology(), ViscousRheology(ν=1000))
        advections = (WENO(), UpwindBiased(order=5))

        ice_thermodynamics = (nothing, SlabSeaIceThermodynamics(grid))

        coriolises = (nothing, FPlane(latitude=45), BetaPlane(latitude=45))
        solvers = (ExplicitSolver(), SplitExplicitSolver())

        for coriolis in coriolises, advection in advections, rheology in rheologies, ice_thermodynamics in ice_thermodynamics, solver in solvers
            dynamics = SeaIceMomentumEquation(grid; coriolis, rheology, solver)

            model = SeaIceModel(grid; dynamics, ice_thermodynamics, advection)
            simulation = Simulation(model, Δt=1.0, stop_iteration=1)

            run!(simulation)

            push!(iterations, model.clock.iteration)
        end
    end

    DistributedComputations.all_reduce!(+, iterations, arch)

    @root jldsave("iterations.jld2"; iterations)
"""

@testset "Sea ice Models" begin
    @info "Testing distributed sea ice models runs..."
    write("distributed_sea_ice_tests.jl", distributed_sea_ice)
    run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project -O0 --check-bounds=yes distributed_sea_ice_tests.jl`)
    rm("distributed_sea_ice_tests.jl")

    @info "Checking that all cores ran all configurations up to 1 iteration..."
    @test isfile("iterations.jld2")
    data = jldopen("iterations.jld2")["iterations"]
    @test all(data .== 4)
    rm("iterations.jld2")
end
