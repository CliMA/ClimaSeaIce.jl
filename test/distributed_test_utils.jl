using JLD2
using MPI
using Oceananigans.DistributedComputations: reconstruct_global_field, reconstruct_global_grid
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics: SeaIceMomentumEquation, ElastoViscoPlasticRheology, SemiImplicitStress

# Run the distributed grid simulation and save down reconstructed results
function run_distributed_sea_ice(arch, filename)
    distributed_grid = RectilinearGrid(arch; 
                                       size = (100, 100, 1), 
                                       x = (-10kilometers, 10kilometers), 
                                       z = (-1, 0), 
                                       halo = (5, 5, 5))

    model = run_distributed_simulation(distributed_grid)

    u = reconstruct_global_field(model.velocities.u)
    v = reconstruct_global_field(model.velocities.v)

    if arch.local_rank == 0
        jldsave(filename; u = Array(interior(u, :, :, 1)),
                          v = Array(interior(v, :, :, 1)))
    end

    MPI.Barrier(MPI.COMM_WORLD)
    MPI.Finalize()

    return nothing
end

# Just a random simulation on a rectilinear grid
function run_distributed_simulation(grid)

    τᵤ = 0.01
    τᵥ = 0.01
    τₒ = SemiImplicitStress()

    # We use an elasto-visco-plastic rheology and WENO seventh order 
    # for advection of h and ℵ
    dynamics = SeaIceMomentumEquation(grid; 
                                      top_momentum_stress = (u=τᵤ, v=τᵥ),
                                      bottom_momentum_stress = τₒ, 
                                      rheology = ElastoViscoPlasticRheology(),
                                      solver = SplitExplicitSolver(grid, substeps=50))

    model = SeaIceModel(; grid, ice_thermodyamics = nothing, dynamics, advection = WENO(order=7))

    for N in 1:100
        time_step!(model, 1minutes)
    end

    return model
end

