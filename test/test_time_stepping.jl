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

        ice_thermodynamics = (nothing, SlabThermodynamics(grid))
        snow_thermodynamics_options = (nothing, SlabSnowThermodynamics(grid))

        coriolises = (nothing, FPlane(latitude=45), BetaPlane(latitude=45))
        solvers = (ExplicitSolver(), SplitExplicitSolver())

        for coriolis in coriolises,
            advection in advections,
            rheology in rheologies,
            ice_thermo in ice_thermodynamics,
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