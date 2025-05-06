function time_step_sea_ice_model_works(grid;
                                       dynamics = nothing,
                                       ice_thermodynamics = nothing,
                                       advection = nothing)

    model = SeaIceModel(grid; dynamics, ice_thermodynamics, advection)
    simulation = Simulation(model, Δt=1.0, stop_iteration=1)

    run!(simulation)

    return model.clock.iteration == 1
end

@testset "Sea ice Models" begin
    @info "Testing sea icemodels..."

    grid_2d = RectilinearGrid(size=(10, 10),     x=(0, 1), y=(0, 1),            topology=(Bounded, Bounded, Flat))
    grid_3d = RectilinearGrid(size=(10, 10, 10), x=(0, 1), y=(0, 1), z=(-1, 0), topology=(Bounded, Bounded, Bounded))

    for grid in (grid_2d, grid_3d)
        rheologies = (ElastoViscoPlasticRheology(), ViscousRheology(ν=1000))
        advections = (WENO(), UpwindBiased(order=5))

        thermodynamics = (nothing, SlabSeaIceThermodynamics(grid))

        coriolises = (nothing, FPlane(latitude=45), BetaPlane(latitude=45))
        solvers = (ExplicitSolver(), SplitExplicitSolver())

        for coriolis in coriolises, advection in advections, rheology in rheologies, ice_thermodynamics in thermodynamics, solver in solvers
            dynamics = SeaIceMomentumEquation(grid; coriolis, rheology, solver)

            @test time_step_sea_ice_model_works(grid;
                                                dynamics,
                                                ice_thermodynamics,
                                                advection)
        end
    end
end