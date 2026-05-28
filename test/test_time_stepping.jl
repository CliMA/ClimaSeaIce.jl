function time_step_sea_ice_model_works(grid;
                                       dynamics = nothing,
                                       ice_thermodynamics = nothing,
                                       snow_thermodynamics = nothing,
                                       advection = nothing)

    model = SeaIceModel(grid; dynamics, ice_thermodynamics, snow_thermodynamics, advection)

    if !isnothing(snow_thermodynamics)
        set!(model, h=1, ℵ=1, hs=0.1)
    end

    simulation = Simulation(model, Δt=1.0, stop_iteration=1)

    run!(simulation)

    return model.clock.iteration == 1
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


function tripolar_evp_diagnostics(fold_topology; Nx = 48, Ny = 24, Δt = 1.0)
    grid = TripolarGrid(CPU(); size = (Nx, Ny, 1), z = (-1, 0), halo = (5, 5, 1), fold_topology = RightCenterFolded)
    solver = SplitExplicitSolver(; substeps = 8)
    dynamics = SeaIceMomentumEquation(grid; rheology = ElastoViscoPlasticRheology(), solver)
    model = SeaIceModel(grid; dynamics, ice_thermodynamics = nothing, advection = nothing)

    set!(model,
         h = 1,
         ℵ = 1,
         u = (x, y, z) -> cosd(x) * exp(-((y - 78) / 7)^2),
         v = (x, y, z) -> sind(2x) * exp(-((y - 78) / 7)^2))

    Oceananigans.BoundaryConditions.fill_halo_regions!(model.velocities.u)
    Oceananigans.BoundaryConditions.fill_halo_regions!(model.velocities.v)

    model_fields = merge(fields(model.dynamics),
                         model.velocities,
                         (; h = model.ice_thickness,
                            ℵ = model.ice_concentration,
                            ρ = model.sea_ice_density))

    rheology = model.dynamics.rheology
    ClimaSeaIce.Rheologies.initialize_rheology!(model, rheology)
    ClimaSeaIce.Rheologies.compute_stresses!(model.dynamics, model_fields, grid, rheology, Δt)
    ClimaSeaIce.Rheologies.finalize_rheology!(model_fields, rheology)

    Sx, Sy, _ = Oceananigans.Grids.worksize(grid)
    kp = ClimaSeaIce.Rheologies.evp_kernel_parameters(grid)

    divu = [ClimaSeaIce.Rheologies.∂ⱼ_σ₁ⱼ(i, j, 1, grid, rheology, model.clock, model_fields) for i in 1:Nx, j in 1:Ny]
    divv = [ClimaSeaIce.Rheologies.∂ⱼ_σ₂ⱼ(i, j, 1, grid, rheology, model.clock, model_fields) for i in 1:Nx, j in 1:Ny]
    combined = max.(abs.(divu), abs.(divv))

    return (; kp,
              Sx,
              Sy,
              row_nm2 = maximum(@view combined[:, Ny-2]),
              row_nm1 = maximum(@view combined[:, Ny-1]),
              row_n = maximum(@view combined[:, Ny]))
end

@testset "RightCenterFolded EVP diagnostics" begin
    diagnostics = rightcenterfolded_evp_diagnostics()

    @test diagnostics.kp == Oceananigans.Utils.KernelParameters(0:diagnostics.Sx+1, 0:diagnostics.Sy+1)

    neighbor_scale = max(diagnostics.row_nm2, diagnostics.row_n, sqrt(eps(typeof(diagnostics.row_nm1))))
    @test diagnostics.row_nm1 <= 5 * neighbor_scale
end
