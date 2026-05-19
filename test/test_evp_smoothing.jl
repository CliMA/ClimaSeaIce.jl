using ClimaSeaIce
using Oceananigans

step_concentration(x, y, z) = x < 0.5 ? 1.0 : 0.05

function initialized_pressure_profile(rheology)
    grid = RectilinearGrid(size=(7, 1), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))
    dynamics = SeaIceMomentumEquation(grid; rheology, solver=SplitExplicitSolver(grid; substeps=1))
    model = SeaIceModel(grid; dynamics, ice_thermodynamics=nothing, advection=WENO())

    set!(model, h=1, ℵ=step_concentration)
    ClimaSeaIce.Rheologies.initialize_rheology!(model, rheology)

    return Array(interior(model.dynamics.auxiliaries.fields.P, :, 1, 1))
end

max_adjacent_pressure_jump(P) = maximum(abs.(diff(P)))

@testset "EVP ice-strength smoothing" begin
    @info "Testing EVP ice-strength smoothing"

    unsmoothed = initialized_pressure_profile(ElastoViscoPlasticRheology(;
        ice_strength_smoothing=0.0,
        ice_strength_smoothing_threshold=1.0,
    ))

    smoothed = initialized_pressure_profile(ElastoViscoPlasticRheology(;
        ice_strength_smoothing=1.0,
        ice_strength_smoothing_threshold=1.0,
    ))

    @test max_adjacent_pressure_jump(smoothed) < max_adjacent_pressure_jump(unsmoothed)
end
