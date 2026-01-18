using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics
using ClimaSeaIce.SeaIceThermodynamics
using Test

using Oceananigans.Fields: @allowscalar
using Oceananigans: prognostic_fields

# Same test as in Oceananigans
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

# Test copied from Oceananigans.jl, adapted for new checkpointing API
function run_checkpointer_tests(true_model, test_model, Δt)
    true_simulation = Simulation(true_model, Δt=Δt, stop_iteration=5)

    checkpointer = Checkpointer(true_model, schedule=IterationInterval(5), overwrite_existing=true)
    true_simulation.output_writers[:checkpointer] = checkpointer

    run!(true_simulation) # for 5 iterations

    checkpointed_model = deepcopy(true_simulation.model)

    true_simulation.stop_iteration = 9
    run!(true_simulation) # for 4 more iterations

    #####
    ##### Test `set!(simulation; checkpoint=filepath)`
    #####

    test_simulation = Simulation(test_model, Δt=Δt, stop_iteration=9)
    test_simulation.output_writers[:checkpointer] =
        Checkpointer(test_model, schedule=IterationInterval(5), overwrite_existing=true)

    set!(test_simulation, checkpoint="checkpoint_iteration5.jld2")

    @test test_simulation.model.clock.iteration == checkpointed_model.clock.iteration
    @test test_simulation.model.clock.time == checkpointed_model.clock.time
    test_model_equality(test_simulation.model, checkpointed_model)

    # This only applies to QuasiAdamsBashforthTimeStepper:
    @test test_simulation.model.clock.last_Δt == checkpointed_model.clock.last_Δt

    #####
    ##### Test pickup from explicit checkpoint path
    #####

    # Recreate test_simulation with fresh model
    test_model2 = deepcopy(test_model)
    test_simulation2 = Simulation(test_model2, Δt=Δt, stop_iteration=9)
    test_simulation2.output_writers[:checkpointer] =
        Checkpointer(test_model2, schedule=IterationInterval(5), overwrite_existing=true)

    # Pickup from explicit checkpoint path
    run!(test_simulation2, pickup="checkpoint_iteration0.jld2")

    @info "Testing model equality when running with pickup=checkpoint_iteration0.jld2."
    @test test_simulation2.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation2.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_simulation2.model, true_simulation.model)

    # Recreate for next test
    test_model3 = deepcopy(test_model)
    test_simulation3 = Simulation(test_model3, Δt=Δt, stop_iteration=9)
    test_simulation3.output_writers[:checkpointer] =
        Checkpointer(test_model3, schedule=IterationInterval(5), overwrite_existing=true)

    run!(test_simulation3, pickup="checkpoint_iteration5.jld2")
    @info "Testing model equality when running with pickup=checkpoint_iteration5.jld2."

    @test test_simulation3.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation3.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_simulation3.model, true_simulation.model)

    #####
    ##### Test `run!(sim, pickup=true)
    #####

    # Recreate for pickup=true test
    test_model4 = deepcopy(test_model)
    test_simulation4 = Simulation(test_model4, Δt=Δt, stop_iteration=9)
    test_simulation4.output_writers[:checkpointer] =
        Checkpointer(test_model4, schedule=IterationInterval(5), overwrite_existing=true)

    run!(test_simulation4, pickup=true)
    @info "    Testing model equality when running with pickup=true."

    @test test_simulation4.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation4.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_simulation4.model, true_simulation.model)

    # Test pickup with iteration number
    test_model5 = deepcopy(test_model)
    test_simulation5 = Simulation(test_model5, Δt=Δt, stop_iteration=9)
    test_simulation5.output_writers[:checkpointer] =
        Checkpointer(test_model5, schedule=IterationInterval(5), overwrite_existing=true)

    run!(test_simulation5, pickup=0)
    @info "    Testing model equality when running with pickup=0."

    @test test_simulation5.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation5.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_simulation5.model, true_simulation.model)

    test_model6 = deepcopy(test_model)
    test_simulation6 = Simulation(test_model6, Δt=Δt, stop_iteration=9)
    test_simulation6.output_writers[:checkpointer] =
        Checkpointer(test_model6, schedule=IterationInterval(5), overwrite_existing=true)

    run!(test_simulation6, pickup=5)
    @info "    Testing model equality when running with pickup=5."

    @test test_simulation6.model.clock.iteration == true_simulation.model.clock.iteration
    @test test_simulation6.model.clock.time == true_simulation.model.clock.time
    test_model_equality(test_simulation6.model, true_simulation.model)

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

@testset "Checkpointing Tests" begin
    test_sea_ice_checkpointer_output(CPU())
end
