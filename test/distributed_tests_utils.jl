using JLD2
using MPI
using Test
using Oceananigans
using Oceananigans.Grids: halo_size, topology
using Oceananigans.Utils: KernelParameters
using Oceananigans.DistributedComputations: reconstruct_global_field, reconstruct_global_grid
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics: SeaIceMomentumEquation, SplitExplicitSolver, ElastoViscoPlasticRheology, SemiImplicitStress

function test_extended_halos()
    arch = Distributed(CPU())
    grid = RectilinearGrid(arch; size=(48, 8, 4), extent=(1, 1, 1),
                        halo=(4, 4, 4),
                        topology=(Periodic, Bounded, Bounded))

    substeps = 8
    solver   = SplitExplicitSolver(grid; substeps)
    dynamics = SeaIceMomentumEquation(grid; solver, rheology=ElastoViscoPlasticRheology(eltype(grid)))
    model    = ClimaSeaIce.SeaIceModel(grid; dynamics)

    Hx, Hy, Hz = halo_size(grid)
    Hx = max(substeps + 2, Hx)
    @test halo_size(model.velocities.u.grid) == (Hx, Hy, Hz)
    @test halo_size(model.dynamics.auxiliaries.fields.P.grid) == (Hx, Hy, Hz)
    # Prognostic tracer-like fields must share the extended halos with the EVP
    # auxiliaries — otherwise the stress kernel reads h/ℵ out-of-bounds.
    @test halo_size(model.ice_thickness.grid)     == (Hx, Hy, Hz)
    @test halo_size(model.ice_concentration.grid) == (Hx, Hy, Hz)

    Nx, Ny, Nz = size(grid)

    kp = model.dynamics.solver.kernel_parameters
    # `split_explicit_kernel_size` for a ConnectedTopology in x is `-H+2:N+H-1`;
    # in y this run is Bounded, so the kernel just iterates `1:N`.
    @test kp == KernelParameters(-Hx+2:Nx+Hx-1, 1:Ny)
end

# Run the distributed grid simulation and save down reconstructed results
function run_distributed_sea_ice(arch, filename)
    distributed_grid = RectilinearGrid(arch; 
                                       size = (100, 100, 1), 
                                       x = (-10kilometers, 10kilometers), 
                                       y = (-10kilometers, 10kilometers), 
                                       z = (-1, 0), 
                                       halo = (5, 5, 5))

    model = run_distributed_simulation(distributed_grid)

    u = reconstruct_global_field(model.velocities.u)
    v = reconstruct_global_field(model.velocities.v)
    h = reconstruct_global_field(model.ice_thickness)
    ℵ = reconstruct_global_field(model.ice_concentration)

    if arch.local_rank == 0
        jldsave(filename; u = Array(interior(u, :, :, 1)),
                          v = Array(interior(v, :, :, 1)),
                          h = Array(interior(h, :, :, 1)),
                          ℵ = Array(interior(ℵ, :, :, 1)))
    end

    MPI.Barrier(MPI.COMM_WORLD)
    MPI.Finalize()

    return nothing
end

# Run the distributed grid simulation and save down reconstructed results
function run_distributed_sea_ice_jl2dwriter(arch, filename)
    distributed_grid = RectilinearGrid(arch; 
                                       size = (100, 100, 1), 
                                       x = (-10kilometers, 10kilometers), 
                                       y = (-10kilometers, 10kilometers), 
                                       z = (-1, 0), 
                                       halo = (5, 5, 5))

    model = run_distributed_jld2_simulation(distributed_grid, filename)

    MPI.Barrier(MPI.COMM_WORLD)
    MPI.Finalize()

    return filename
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

    model = SeaIceModel(grid; dynamics, advection = WENO(order=7))

    for N in 1:100
        time_step!(model, 1minutes)
    end

    return model
end

# Just a random simulation on a rectilinear grid
function run_distributed_jld2_simulation(grid, filename)

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

    model = SeaIceModel(grid; dynamics, advection = WENO(order=7))
    simulation = Simulation(model, Δt = 10, stop_iteration = 20)

    fields = Oceananigans.prognostic_fields(model)
    simulation.output_writers[:fields] = JLD2Writer(model, fields;
                                                    schedule = IterationInterval(1),
                                                    filename,
                                                    overwrite_existing = true)

    run!(simulation)
    
    return model
end
