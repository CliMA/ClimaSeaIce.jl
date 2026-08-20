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

@testset "Ocean surface tilt" begin
    @info "Testing the ocean surface tilt for explicit and split-explicit solvers..."

    L = 10_000
    grid = RectilinearGrid(size = (8, 8, 1), x = (0, L), y = (0, L),
                           z = (-1, 0), halo = (4, 4, 4), topology = (Periodic, Periodic, Bounded))

    g = Oceananigans.defaults.gravitational_acceleration
    Δt = 60

    η = Field{Center, Center, Nothing}(grid)
    set!(η, (x, y) -> 0.1 * sin(2π * x / L) * sin(2π * y / L))

    ∂xη = Field(∂x(η))
    ∂yη = Field(∂y(η))
    compute!(∂xη)
    compute!(∂yη)

    for solver in (ExplicitSolver(), SplitExplicitSolver(grid; substeps=10))
        ocean_surface_tilt = OceanSurfaceTilt(; η, gravitational_acceleration = g)
        dynamics = SeaIceMomentumEquation(grid; ocean_surface_tilt, rheology = ViscousRheology(ν=0), solver)
        model = SeaIceModel(grid; dynamics)
        set!(model, h = 1, ℵ = 1, u = 0, v = 0)

        time_step!(model, Δt)

        u = interior(model.velocities.u, :, :, 1)
        v = interior(model.velocities.v, :, :, 1)

        # Starting from rest, the ice accelerates down the ocean surface slope
        @test u ≈ - g * Δt * interior(∂xη, :, :, 1)
        @test v ≈ - g * Δt * interior(∂yη, :, :, 1)
    end
end
