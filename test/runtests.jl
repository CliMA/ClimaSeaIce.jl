using Test

using Oceananigans
using ClimaSeaIce
using ClimaSeaIce.EnthalpyMethodSeaIceModels: EnthalpyMethodSeaIceModel, MolecularDiffusivity

TEST_GROUP = get(ENV, "TEST_GROUP", "all")

if TEST_GROUP == "all" || TEST_GROUP == "enthalpy"
    κ = 1e-5
    grid = RectilinearGrid(size=3, z=(-3, 0), topology=(Flat, Flat, Bounded))
    closure = MolecularDiffusivity(grid, κ_ice=κ, κ_water=κ)
    model = EnthalpyMethodSeaIceModel(; grid, closure)

    @test model isa EnthalpyMethodSeaIceModel

    # Test that it runs
    simulation = Simulation(model; Δt = 0.1 / κ, stop_iteration=3)

    @test begin
        try
            run!(simulation)
            true
        catch
            false
        end
    end
end

if TEST_GROUP == "all" || TEST_GROUP == "advection"
    include("test_sea_ice_advection.jl")
end

if TEST_GROUP == "all" || TEST_GROUP == "timestepping"
    include("test_time_stepping.jl")
end

if TEST_GROUP == "all" || TEST_GROUP == "checkpointing"
    include("test_checkpointing.jl")
end

if TEST_GROUP == "all" || TEST_GROUP == "distributed"
    include("test_distributed_sea_ice.jl")
end
