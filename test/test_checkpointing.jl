using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics
using ClimaSeaIce.SeaIceThermodynamics
using Test

using Oceananigans.Fields: @allowscalar
using Oceananigans: prognostic_fields

# # Same test as in Oceananigans
function test_model_equality(test_model, true_model)
    @allowscalar begin
        test_model_fields = prognostic_fields(test_model)
        true_model_fields = prognostic_fields(true_model)
        field_names = keys(test_model_fields)

        for name in field_names
            @test all(test_model_fields[name].data .≈ true_model_fields[name].data)
        end
    end

    return nothing
end

# function initialization_test_simulation(arch, stop_time, Δt=1, δt=2)
#     grid = RectilinearGrid(arch, size=(), topology=(Flat, Flat, Flat))
#     model = SeaIceModel(grid)
#     simulation = Simulation(model; Δt, stop_time)

#     progress_message(sim) = @info string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
#     simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(δt))

#     checkpointer = Checkpointer(model,
#                                 schedule = TimeInterval(stop_time),
#                                 prefix = "initialization_test",
#                                 cleanup = false)

#     simulation.output_writers[:checkpointer] = checkpointer

#     return simulation
# end

# Same test as in Oceananigans
function run_checkpointer_tests(true_model, test_model, Δt)
    true_simulation = Simulation(true_model, Δt=Δt, stop_iteration=5)

    checkpointer = Checkpointer(true_model, schedule=IterationInterval(5), overwrite_existing=true)
    push!(true_simulation.output_writers, checkpointer)

    run!(true_simulation) # for 5 iterations

    checkpointed_model = deepcopy(true_simulation.model)

    true_simulation.stop_iteration = 9
    run!(true_simulation) # for 4 more iterations

    #####
    ##### Test `set!(model, checkpoint_file)`
    #####

    set!(test_model, "checkpoint_iteration5.jld2")

    @test test_model.clock.iteration == checkpointed_model.clock.iteration
    @test test_model.clock.time == checkpointed_model.clock.time
    test_model_equality(test_model, checkpointed_model)

    # This only applies to QuasiAdamsBashforthTimeStepper:
    @test test_model.clock.last_Δt == checkpointed_model.clock.last_Δt

    #####
    ##### Test pickup from explicit checkpoint path
    #####

    test_simulation = Simulation(test_model, Δt=Δt, stop_iteration=9)

    # Pickup from explicit checkpoint path
    run!(test_simulation, pickup="checkpoint_iteration0.jld2")

    @info "Testing model equality when running with pickup=checkpoint_iteration0.jld2."
    @test test_simulation.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_model, true_model)

    run!(test_simulation, pickup="checkpoint_iteration5.jld2")
    @info "Testing model equality when running with pickup=checkpoint_iteration5.jld2."

    @test test_simulation.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_model, true_model)

    #####
    ##### Test `run!(sim, pickup=true)
    #####

    # Pickup using existing checkpointer
    test_simulation.output_writers[:checkpointer] =
        Checkpointer(test_model, schedule=IterationInterval(5), overwrite_existing=true)

    run!(test_simulation, pickup=true)
    @info "    Testing model equality when running with pickup=true."

    @test test_simulation.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_model, true_model)

    run!(test_simulation, pickup=0)
    @info "    Testing model equality when running with pickup=0."

    @test test_simulation.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_model, true_model)

    run!(test_simulation, pickup=5)
    @info "    Testing model equality when running with pickup=5."

    @test test_simulation.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_model, true_model)

    rm("checkpoint_iteration0.jld2", force=true)
    rm("checkpoint_iteration5.jld2", force=true)

    return nothing
end

function test_sea_ice_checkpointer_output(arch)
    # Create and run "true model"
    Nx, Ny = 16, 16
    Lx, Ly = 100, 100
    Δt = 1

    grid = RectilinearGrid(arch, size=(Nx, Ny), x=(0, Lx), y=(0, Ly), topology=(Bounded, Bounded, Flat))
    for ice_thermodynamics in (nothing, SlabSeaIceThermodynamics(grid))
        for dynamics in (nothing, SeaIceMomentumEquation(grid))
            
            true_model = SeaIceModel(grid; dynamics, ice_thermodynamics)
            test_model = deepcopy(true_model)

            for field in merge(true_model.velocities, 
                            (; h = true_model.ice_concentration, 
                               ℵ = true_model.ice_thickness))

                set!(field, (x, y) -> rand() * 1e-5)
            end

            run_checkpointer_tests(true_model, test_model, Δt)
        end
    end
end