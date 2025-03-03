i

function time_step_sea_ice_model_works(grid;
                                       dynamics = nothing,
                                       ice_thermodynamics = nothing,
                                       advection = nothing)

    model = SeaIceModel(; grid, dynamics, ice_thermodynamics, advection)
    simulation = Simulation(model, Δt=1.0, stop_iteration=1)

    run!(simulation)

    return model.clock.iteration == 1
end

@testset "Sea ice Models" begin
    @info "Testing sea icemodels..."

    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))

    rheologies = (ElastoViscoPlasticRheology(), ViscousRheology(ν=1000))
    advections = (WENO(), UpwindBiased(order=5))

    ice_thermodynamics = (nothing, SlabSeaIceThermodynamics())

    coriolises = (nothing, BetaPlane())
    solvers = (ExplicitSolver(), SpliExplicitSolver())

    for coriolis in coriolises, advection in advections, rheology in rheologies, thermodynamics in ice_thermodynamics, solver in solvers
        dynamics = SeaIceMomentumEquation(grid; coriolis, rheology, solver)

        @test time_step_sea_ice_model_works(grid;
                                            dynamics,
                                            ice_thermodynamics=thermodynamics,
                                            advection=advection)
    end

end