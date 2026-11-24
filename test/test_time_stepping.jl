using Oceananigans.Fields: ConstantField

function time_step_sea_ice_model_works(grid;
                                       dynamics = nothing,
                                       ice_thermodynamics = nothing,
                                       snow_thermodynamics = nothing,
                                       advection = nothing)

    model = SeaIceModel(grid; dynamics, ice_thermodynamics, snow_thermodynamics, advection)

    if !isnothing(snow_thermodynamics)
        set!(model, h=1, ℵ=1, hs=0.1)
    end

    simulation = Simulation(model, Δt=1.1, stop_iteration=1, verbose=false)

    run!(simulation)

    return model.clock.iteration == 1 && model.clock.time == 1.1
end

@testset "Sea ice Models" begin
    @info "Testing sea ice models..."

    grid_2d = RectilinearGrid(size=(10, 10),     x=(0, 1), y=(0, 1),            topology=(Bounded, Bounded, Flat))
    grid_3d = RectilinearGrid(size=(10, 10, 10), x=(0, 1), y=(0, 1), z=(-1, 0), topology=(Bounded, Bounded, Bounded))

    for grid in (grid_2d, grid_3d)
        rheologies = (ElastoViscoPlasticRheology(), ViscousRheology(ν=1000))
        advections = (WENO(), UpwindBiased(order=5))

        ice_thermodynamics_options = (nothing, SlabThermodynamics(grid))
        snow_thermodynamics_options = (nothing, snow_slab_thermodynamics(grid))

        coriolises = (nothing, FPlane(latitude=45), BetaPlane(latitude=45))
        solvers = (ExplicitSolver(), SplitExplicitSolver())

        for coriolis in coriolises,
            advection in advections,
            rheology in rheologies,
            ice_thermo in ice_thermodynamics_options,
            snow_thermo in snow_thermodynamics_options,
            solver in solvers

            dynamics = SeaIceMomentumEquation(grid; coriolis, rheology, solver)

            @test time_step_sea_ice_model_works(grid;
                                                dynamics,
                                                ice_thermodynamics = ice_thermo,
                                                snow_thermodynamics = snow_thermo,
                                                advection)
        end
    end
end

@testset "Sea ice dynamics without a rheology" begin
    @info "Testing sea ice dynamics without a rheology..."

    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))

    for solver in (ExplicitSolver(), SplitExplicitSolver())
        dynamics = SeaIceMomentumEquation(grid; rheology = nothing, solver)
        @test time_step_sea_ice_model_works(grid; dynamics)
    end
end

@testset "SemiImplicitStress ocean drag" begin
    @info "Testing SemiImplicitStress ocean drag for explicit and split-explicit solvers..."

    grid = RectilinearGrid(size = (8, 8, 1), x = (0, 10_000), y = (0, 10_000),
                           z = (-1, 0), halo = (4, 4, 4), topology = (Periodic, Periodic, Bounded))

    uₒ = 0.1 # ocean x-velocity dragging the ice from rest

    for solver in (ExplicitSolver(), SplitExplicitSolver(grid; substeps=10))
        τₒ = SemiImplicitStress(uₑ = ConstantField(uₒ))
        dynamics = SeaIceMomentumEquation(grid; bottom_momentum_stress = τₒ,
                                          rheology = ElastoViscoPlasticRheology(), solver)
        model = SeaIceModel(grid; dynamics)
        set!(model, h = 1, ℵ = 1, u = 0, v = 0)

        for _ in 1:20
            time_step!(model, 60)
        end

        u = interior(model.velocities.u)
        @test all(isfinite, u)
        @test maximum(u) > 0   # the drag accelerated the ice toward the moving ocean
        @test maximum(u) ≤ uₒ  # but never overshoots the ocean velocity
    end
end

@testset "SemiImplicitStress computes its external velocities" begin
    @info "Testing that computed external velocities follow their source..."

    grid = RectilinearGrid(size = (8, 8, 1), x = (0, 10_000), y = (0, 10_000),
                           z = (-1, 0), halo = (4, 4, 4), topology = (Periodic, Periodic, Bounded))

    for solver in (ExplicitSolver(), SplitExplicitSolver(grid; substeps=10))
        uₒ = Field{Face, Center, Nothing}(grid)
        τₒ  = SemiImplicitStress(uₑ = Field(2 * uₒ))
        dynamics = SeaIceMomentumEquation(grid; bottom_momentum_stress = τₒ,
                                          rheology = ElastoViscoPlasticRheology(), solver)
        model = SeaIceModel(grid; dynamics)
        set!(model, h = 1, ℵ = 1, u = 0, v = 0)
        set!(uₒ, 0.05)

        time_step!(model, 60)

        @test maximum(interior(model.velocities.u)) > 0
    end
end

function tripolar_evp_diagnostics(fold_topology; Nx = 48, Ny = 24, Δt = 1.0)
    grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), z = (-1, 0), fold_topology = fold_topology)
    Hx, Hy, _ = Oceananigans.Grids.halo_size(grid)
    solver = SplitExplicitSolver(; substeps = 8)
    dynamics = SeaIceMomentumEquation(grid; rheology = ElastoViscoPlasticRheology(), solver)
    model = SeaIceModel(grid; dynamics, ice_thermodynamics = nothing, advection = nothing)

    set!(model,
     h = 1,
     ℵ = 1,
     u = (x, y) -> cosd(x) * exp(-((y - 78) / 7)^2),
     v = (x, y) -> sind(2x) * exp(-((y - 78) / 7)^2))
  
    Oceananigans.BoundaryConditions.fill_halo_regions!(model.velocities.u)
    Oceananigans.BoundaryConditions.fill_halo_regions!(model.velocities.v)

    model_fields = merge(fields(model.dynamics),
                         model.velocities,
                         (; h = model.ice_thickness,
                            ℵ = model.ice_concentration,
                            ρ = model.sea_ice_density))

    rheology = model.dynamics.rheology
    ClimaSeaIce.Rheologies.initialize_rheology!(model, rheology)
    u_immersed_bc = model.velocities.u.boundary_conditions.immersed
    v_immersed_bc = model.velocities.v.boundary_conditions.immersed
    ClimaSeaIce.Rheologies.compute_stresses!(model.dynamics, model_fields, grid, rheology, Δt, u_immersed_bc, v_immersed_bc)
    ClimaSeaIce.Rheologies.finalize_rheology!(model_fields, rheology)

    Sx, Sy, _ = Oceananigans.Utils.worksize(grid)
    kp = ClimaSeaIce.Rheologies.evp_kernel_parameters(grid)

    divu = [ClimaSeaIce.Rheologies.∂ⱼ_σ₁ⱼ(i, j, 1, grid, rheology, model.clock, model_fields, u_immersed_bc) for i in 1:Nx, j in 1:Ny]
    divv = [ClimaSeaIce.Rheologies.∂ⱼ_σ₂ⱼ(i, j, 1, grid, rheology, model.clock, model_fields, v_immersed_bc) for i in 1:Nx, j in 1:Ny]
    combined = max.(abs.(divu), abs.(divv))

    return (; kp,
              Sx,
              Sy,
              Hx,
              Hy,
              row_nm2 = maximum(@view combined[:, Ny-2]),
              row_nm1 = maximum(@view combined[:, Ny-1]),
              row_n = maximum(@view combined[:, Ny]))
end

@testset "RightFaceFolded EVP diagnostics" begin
    diagnostics = tripolar_evp_diagnostics(Oceananigans.Grids.RightFaceFolded)

    @test diagnostics.kp == Oceananigans.Utils.KernelParameters(-diagnostics.Hx+2:diagnostics.Sx+diagnostics.Hx-1, -diagnostics.Hy+2:diagnostics.Sy+diagnostics.Hy-1)

    neighbor_scale = max(diagnostics.row_nm2, diagnostics.row_n, sqrt(eps(typeof(diagnostics.row_nm1))))
    @test diagnostics.row_nm1 <= 5 * neighbor_scale
end
