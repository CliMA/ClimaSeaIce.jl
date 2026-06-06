using JLD2
using MPI
using Test
using Oceananigans
using Oceananigans.Grids: halo_size, topology
using Oceananigans.Utils: KernelParameters
using Oceananigans.DistributedComputations: @root, ranks, reconstruct_global_field, reconstruct_global_grid
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
    Hx = max(2substeps + 3, Hx)
    @test halo_size(model.velocities.u.grid) == (Hx, Hy, Hz)
    @test halo_size(model.dynamics.auxiliaries.fields.P.grid) == (Hx, Hy, Hz)
    @test halo_size(model.ice_thickness.grid)     == (Hx, Hy, Hz)
    @test halo_size(model.ice_concentration.grid) == (Hx, Hy, Hz)

    Nx, Ny, Nz = size(grid)

    kp = model.dynamics.solver.kernel_parameters
    @test kp == KernelParameters(-Hx+2:Nx+Hx-1, 1:Ny)
end

function test_distributed_simulations()
    grid = RectilinearGrid(CPU();
                           size = (100, 100, 1),
                           x = (-10kilometers, 10kilometers),
                           y = (-10kilometers, 10kilometers),
                           z = (-1, 0),
                           halo = (5, 5, 5))

    models = run_distributed_simulation(grid)

    # Retrieve Serial quantities
    us, vs = models.velocities
    hs = models.ice_thickness
    ℵs = models.ice_concentration

    us = interior(us, :, :, 1)
    vs = interior(vs, :, :, 1)
    hs = interior(hs, :, :, 1)
    ℵs = interior(ℵs, :, :, 1)

    archs = [Distributed(CPU(), partition = Partition(1, 4))
             Distributed(CPU(), partition = Partition(4, 1))
             Distributed(CPU(), partition = Partition(2, 2))]

    for arch in archs
        @root @info "Testing $(ranks(arch)) distributed simulation"
        grid = RectilinearGrid(arch;
                               size = (100, 100, 1),
                               x = (-10kilometers, 10kilometers),
                               y = (-10kilometers, 10kilometers),
                               z = (-1, 0),
                               halo = (5, 5, 5))

        modelp = run_distributed_simulation(grid)
        up, vp = models.velocities
        hp = models.ice_thickness
        ℵp = models.ice_concentration

        up = interior(reconstruct_global_field(up), :, :, 1)
        vp = interior(reconstruct_global_field(vp), :, :, 1)
        hp = interior(reconstruct_global_field(hp), :, :, 1)
        ℵp = interior(reconstruct_global_field(ℵp), :, :, 1)

        @test all(us .≈ up)
        @test all(vs .≈ vp)
        @test all(hs .≈ hp)
        @test all(ℵs .≈ ℵp)
    end
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
                                      coriolis = FPlane(f = 1e-4),
                                      solver = SplitExplicitSolver(grid, substeps=10))

    model = SeaIceModel(grid; dynamics, advection = WENO(order=7))

    # Initialize with non-trivial values
    h₀(x, y) = 0.3 + 0.005 * (sin(60 * x / 20kilometers) + sin(30 * y / 20kilometers))

    set!(model, h = h₀)
    set!(model, ℵ = 1)
    set!(model, u = 0.1)

    for N in 1:20
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
                                      solver = SplitExplicitSolver(grid, substeps=10))

    model = SeaIceModel(grid; dynamics, advection = WENO(order=7))
    simulation = Simulation(model, Δt = 10, stop_iteration = 20, verbose = false)

    fields = Oceananigans.prognostic_fields(model)
    simulation.output_writers[:fields] = JLD2Writer(model, fields;
                                                    schedule = IterationInterval(1),
                                                    filename,
                                                    with_halos = false,
                                                    overwrite_existing = true)

    run!(simulation)

    return model
end
