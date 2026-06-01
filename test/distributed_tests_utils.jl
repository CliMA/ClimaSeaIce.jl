using JLD2
using MPI
using Test
using Oceananigans
using Oceananigans.Grids: halo_size, topology
using Oceananigans.Utils: KernelParameters
using Oceananigans.DistributedComputations: @root, ranks, reconstruct_global_field, reconstruct_global_grid
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.OrthogonalSphericalShellGrids: TripolarGrid
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics: SeaIceMomentumEquation, SplitExplicitSolver, ElastoViscoPlasticRheology,
                                  SemiImplicitStress, StressBalanceFreeDrift, maybe_extended_grid

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
        up, vp = modelp.velocities
        hp = modelp.ice_thickness
        ℵp = modelp.ice_concentration

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

    # We use an elasto-visco-plastic rheology and WENO seventh order for advection of h and ℵ
    dynamics = SeaIceMomentumEquation(grid;
                                      top_momentum_stress = (u=τᵤ, v=τᵥ),
                                      bottom_momentum_stress = τₒ,
                                      rheology = ElastoViscoPlasticRheology(),
                                      solver = SplitExplicitSolver(grid, substeps=10))

    model = SeaIceModel(grid; dynamics, advection = WENO(order=7))
    set!(model, h = 1, ℵ = 0.5)

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

# Mask the two tripolar singularities (within `radius` degrees of each north pole) and the
# southernmost row. Their metric is degenerate, and unmasked the EVP rheology divides by it →
# NaN. A realistic run masks them via bathymetry; here we do it analytically.
function analytical_immersed_tripolar_grid(underlying_grid::TripolarGrid; radius = 5)
    λp = underlying_grid.conformal_mapping.first_pole_longitude
    φp = underlying_grid.conformal_mapping.north_poles_latitude
    φm = underlying_grid.conformal_mapping.southernmost_latitude
    Lz = underlying_grid.Lz

    bottom_height(λ, φ) = ((abs(λ - λp) < radius)       & (abs(φp - φ) < radius)) |
                          ((abs(λ - λp - 180) < radius) & (abs(φp - φ) < radius)) |
                          ((abs(λ - λp - 360) < radius) & (abs(φp - φ) < radius)) | (φ < φm + radius) ? 0 : -Lz

    return ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))
end

function run_tripolar_sea_ice_simulation(grid)
    solver = SplitExplicitSolver(grid; substeps=10)

    # Atmosphere stress is owned here and written externally each step. Build it on the (halo-extended)
    # velocity grid so the substep solver aliases it — no copy, so the live reference is preserved.
    velocity_grid = maybe_extended_grid(solver, grid)
    τua = Field{Face, Center, Nothing}(velocity_grid; boundary_conditions = ClimaSeaIce.default_sea_ice_boundary_conditions(velocity_grid, :u))
    τva = Field{Center, Face, Nothing}(velocity_grid; boundary_conditions = ClimaSeaIce.default_sea_ice_boundary_conditions(velocity_grid, :v))

    # Ocean velocities model an ocean-owned source on the base grid (the SemiImplicitStress keeps the
    # source and refreshes an extended copy each step).
    uₒ  = Field{Face, Center, Nothing}(grid)
    vₒ  = Field{Center, Face, Nothing}(grid)

    set!(τua, (λ, φ) ->  0.05 * sind(λ))
    set!(τva, (λ, φ) ->  0.05 * cosd(λ))
    set!(uₒ,  (λ, φ) ->  0.02 * cosd(λ))
    set!(vₒ,  (λ, φ) -> -0.02 * sind(λ))
    foreach(fill_halo_regions!, (τua, τva, uₒ, vₒ))

    τₒ = SemiImplicitStress(uₑ=uₒ, vₑ=vₒ)

    dynamics = SeaIceMomentumEquation(velocity_grid;
                                      top_momentum_stress = (u=τua, v=τva),
                                      bottom_momentum_stress = τₒ,
                                      rheology = ElastoViscoPlasticRheology(),
                                      coriolis = FPlane(f = 1e-4),
                                      free_drift = StressBalanceFreeDrift((u=τua, v=τva), τₒ),
                                      solver)

    # No tracer advection, enable when https://github.com/CliMA/Oceananigans.jl/pull/5564 is merged
    model = SeaIceModel(grid; dynamics, advection = nothing)

    # The atmosphere stress must not have been copied — the external writer's reference is `τua`.
    @assert model.dynamics.external_momentum_stresses.top.u === τua "atmosphere stress was copied; live reference lost"

    set!(model, h = (λ, φ) -> 0.3 + 0.05 * sind(2λ) * cosd(φ))
    set!(model, ℵ = 1)
    set!(model, u = 0.05)

    for _ in 1:8
        time_step!(model, 1minutes)
    end

    return model
end

# Equality of a distributed (halo-extended) tripolar sea-ice run with the serial (non-extended) one.
function test_distributed_tripolar_simulations()
    serial_grid = analytical_immersed_tripolar_grid(TripolarGrid(CPU(); size = (60, 120, 1), southernmost_latitude = 60, halo = (4, 4, 4), z = (-1, 0)))
    serial = run_tripolar_sea_ice_simulation(serial_grid)

    us = interior(serial.velocities.u, :, :, 1)
    vs = interior(serial.velocities.v, :, :, 1)
    hs = interior(serial.ice_thickness, :, :, 1)
    ℵs = interior(serial.ice_concentration, :, :, 1)

    # Pure y-partition: the fold lives on the northernmost rank and the wide halo crosses ranks in y
    # — the same decomposition as the OMIP run.
    arch = Distributed(CPU(), partition = Partition(1, 4))
    @root @info "Testing $(ranks(arch)) distributed tripolar sea ice simulation equality"

    dist_grid = analytical_immersed_tripolar_grid(TripolarGrid(arch; size = (60, 120, 1), southernmost_latitude = 60, halo = (4, 4, 4), z = (-1, 0)))
    distributed = run_tripolar_sea_ice_simulation(dist_grid)

    up = interior(reconstruct_global_field(distributed.velocities.u), :, :, 1)
    vp = interior(reconstruct_global_field(distributed.velocities.v), :, :, 1)
    hp = interior(reconstruct_global_field(distributed.ice_thickness), :, :, 1)
    ℵp = interior(reconstruct_global_field(distributed.ice_concentration), :, :, 1)

    @test all(us .≈ up)
    @test all(vs .≈ vp)
    @test all(hs .≈ hp)
    @test all(ℵs .≈ ℵp)

    return nothing
end
