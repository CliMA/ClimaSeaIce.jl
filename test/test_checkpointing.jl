using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics
using ClimaSeaIce.SeaIceThermodynamics
using Test

using Oceananigans
using Oceananigans.Fields: @allowscalar, ConstantField
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
    for ice_thermodynamics in (nothing, SlabThermodynamics(grid))
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

@testset "Checkpointing with snow" begin
    Nx, Ny = 16, 16
    Lx, Ly = 100, 100
    Δt = 1

    grid = RectilinearGrid(CPU(), size=(Nx, Ny), x=(0, Lx), y=(0, Ly), topology=(Bounded, Bounded, Flat))
    snow_thermo = snow_slab_thermodynamics(grid)

    true_model = SeaIceModel(grid; ice_thermodynamics=SlabThermodynamics(grid), snow_thermodynamics=snow_thermo)
    test_model = deepcopy(true_model)

    for field in merge(true_model.velocities,
                       (; h  = true_model.ice_thickness,
                          ℵ  = true_model.ice_concentration,
                          hs = true_model.snow_thickness))
        set!(field, (x, y) -> rand() * 1e-5)
    end

    run_checkpointer_tests(true_model, test_model, Δt)
end

@testset "Checkpointing SemiImplicitStress caches" begin
    # Regression test: `SemiImplicitStress` (used for ice-ocean and
    # ice-atmosphere drag) caches per-cell implicit-stress coefficients in
    # `τᵢᵤ`/`τᵢᵥ` Fields. They are written by `compute_implicit_stress_coefficients!`
    # inside the EVP substep loop and read by `_compute_sea_ice_ocean_stress!`
    # (via `implicit_τx_coefficient`) every time the coupled `update_state!`
    # runs — including the call triggered by `initialize_simulation!` on a
    # restored simulation. If `prognostic_state(::SeaIceMomentumEquation)`
    # doesn't walk `external_momentum_stresses`, these caches do not
    # round-trip and the first post-restore flux computation reads zeros,
    # producing ocean-velocity drift in a coupled simulation.

    Δt = 1
    grid = RectilinearGrid(CPU(), size=(16, 16), x=(0, 100), y=(0, 100),
                           topology=(Bounded, Bounded, Flat))

    # SeaIceModel with `SemiImplicitStress` on both ice surfaces — exercises
    # both top and bottom branches of `external_momentum_stresses`.
    function make_model(grid)
        top    = SemiImplicitStress(; uₑ=ConstantField(10), vₑ=ConstantField(5),
                                      ρₑ=1.225,  Cᴰ=1.5e-3)
        bottom = SemiImplicitStress(; uₑ=ConstantField(0),  vₑ=ConstantField(0),
                                      ρₑ=1026.0, Cᴰ=5.5e-3)
        dynamics = SeaIceMomentumEquation(grid; top_momentum_stress    = top,
                                                bottom_momentum_stress = bottom)
        SeaIceModel(grid; dynamics, ice_thermodynamics=SlabThermodynamics(grid))
    end

    function set_random_initial_conditions!(model)
        for field in merge(model.velocities,
                           (h = model.ice_thickness, ℵ = model.ice_concentration))
            set!(field, (x, y) -> rand() * 1e-5)
        end
    end

    # Per-cache snapshot helper for the τᵢᵤ/τᵢᵥ Fields.
    cache(model, side, name) =
        Array(parent(getproperty(getproperty(model.dynamics.external_momentum_stresses, side), name)))

    # Outer round-trip via the existing helper — `prognostic_fields` equality,
    # for `SeaIceMomentumEquation`-with-`SemiImplicitStress` configurations.
    true_model = make_model(grid)
    test_model = deepcopy(true_model)
    set_random_initial_conditions!(true_model)
    run_checkpointer_tests(true_model, test_model, Δt)

    # Focused cache round-trip. `test_model_equality` only covers
    # `prognostic_fields(::SeaIceModel)`, which does not include
    # `external_momentum_stresses` — so the τᵢᵤ/τᵢᵥ round-trip must be checked
    # directly here.
    @testset "τᵢᵤ/τᵢᵥ cache round-trip" begin
        prefix = "semi_implicit_stress_checkpoint"

        # Save: run the EVP loop long enough to populate the caches.
        saved = make_model(grid)
        set_random_initial_conditions!(saved)
        sim_save = Simulation(saved; Δt, stop_iteration=5)
        sim_save.output_writers[:checkpointer] =
            Checkpointer(saved; schedule=IterationInterval(5),
                                overwrite_existing=true, prefix)
        run!(sim_save)

        saved_caches = (top    = (τᵢᵤ = cache(saved, :top, :τᵢᵤ),
                                  τᵢᵥ = cache(saved, :top, :τᵢᵥ)),
                        bottom = (τᵢᵤ = cache(saved, :bottom, :τᵢᵤ),
                                  τᵢᵥ = cache(saved, :bottom, :τᵢᵥ)))

        # Sanity: the substep loop actually populated the caches — otherwise
        # the round-trip is vacuously satisfied by both sides being zero.
        @test maximum(abs.(saved_caches.top.τᵢᵤ))    > 0
        @test maximum(abs.(saved_caches.bottom.τᵢᵤ)) > 0

        # Restore: fresh model, pull state from disk.
        restored = make_model(grid)
        sim_restore = Simulation(restored; Δt, stop_iteration=5)
        sim_restore.output_writers[:checkpointer] =
            Checkpointer(restored; schedule=IterationInterval(5),
                                   overwrite_existing=true, prefix)
        set!(sim_restore, checkpoint="$(prefix)_iteration5.jld2")

        # Round-trip equality must be bit-exact: every τᵢᵤ/τᵢᵥ Field is a
        # deterministic function of the saved prognostic state.
        @test cache(restored, :top,    :τᵢᵤ) == saved_caches.top.τᵢᵤ
        @test cache(restored, :top,    :τᵢᵥ) == saved_caches.top.τᵢᵥ
        @test cache(restored, :bottom, :τᵢᵤ) == saved_caches.bottom.τᵢᵤ
        @test cache(restored, :bottom, :τᵢᵥ) == saved_caches.bottom.τᵢᵥ

        for f in filter(p -> startswith(p, prefix), readdir())
            rm(f; force=true)
        end
    end
end
