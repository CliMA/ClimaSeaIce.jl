using Test
using Oceananigans
using Oceananigans.Fields: ConstantField
using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics: SplitExplicitSolver, SemiImplicitStress, SeaIceMomentumEquation
using ClimaSeaIce.Rheologies: ElastoViscoPlasticRheology

# An ice patch surrounded by open water, dragged around by an ocean eddy, so that
# the ice edge (and therefore the set of skipped cells) changes between time steps
function ice_patch_model(grid; skip_ice_free_cells, rheology = ElastoViscoPlasticRheology(), ice = :patch)
    uₒ = XFaceField(grid)
    vₒ = YFaceField(grid)
    set!(uₒ, (x, y, z...) -> - 0.2 * sin(π * y / 100e3))
    set!(vₒ, (x, y, z...) ->   0.2 * sin(π * x / 100e3))

    τₒ = SemiImplicitStress(uₑ = uₒ, vₑ = vₒ)
    solver = SplitExplicitSolver(grid; substeps = 50, skip_ice_free_cells)
    dynamics = SeaIceMomentumEquation(grid; coriolis = FPlane(latitude = 70),
                                      bottom_momentum_stress = τₒ, rheology, solver)

    model = SeaIceModel(grid; dynamics, advection = WENO(order = 5))

    patch(x, y) = (30e3 < x < 60e3) & (40e3 < y < 80e3)
    h = ice == :patch ? (x, y, z...) -> patch(x, y) * (1 + x / 100e3) :
        ice == :none  ? 0 : 1
    ℵ = ice == :patch ? (x, y, z...) -> patch(x, y) * 0.9 :
        ice == :none  ? 0 : 0.9

    set!(model, h = h, ℵ = ℵ)
    return model
end

function compare_skipping(grid; Nt = 10, perturb_velocities = false, kw...)
    full = ice_patch_model(grid; skip_ice_free_cells = false, kw...)
    skip = ice_patch_model(grid; skip_ice_free_cells = true,  kw...)

    for n in 1:Nt
        # Velocities set over open water must be zeroed by the first substep
        if perturb_velocities && n == Nt ÷ 2
            for model in (full, skip)
                set!(model, u = 0.05, v = -0.05)
            end
        end

        time_step!(full, 600)
        time_step!(skip, 600)
    end

    fields_to_compare = (:u, :v, :h, :ℵ)
    same = Dict{Symbol, Bool}()
    for name in fields_to_compare
        same[name] = parent(fields(full)[name]) == parent(fields(skip)[name])
    end

    aux_full = full.dynamics.auxiliaries.fields
    aux_skip = skip.dynamics.auxiliaries.fields
    for name in keys(aux_full)
        name == :Δ && continue # diagnostic, left stale in ice-free cells
        same[name] = parent(aux_full[name]) == parent(aux_skip[name])
    end

    return full, skip, same
end

@testset "Skipping ice-free cells in the split-explicit solver" begin
    periodic = RectilinearGrid(size = (40, 40), x = (0, 100e3), y = (0, 100e3), halo = (7, 7),
                               topology = (Periodic, Periodic, Flat))

    bounded = RectilinearGrid(size = (40, 40), x = (0, 100e3), y = (0, 100e3), halo = (7, 7),
                              topology = (Bounded, Bounded, Flat))

    # An island in the middle of the patch
    underlying = RectilinearGrid(size = (40, 40, 1), x = (0, 100e3), y = (0, 100e3), z = (-10, 0),
                                 halo = (7, 7, 7), topology = (Bounded, Bounded, Bounded))
    island(x, y) = ifelse((45e3 < x < 55e3) & (55e3 < y < 65e3), 0, -10)
    immersed = ImmersedBoundaryGrid(underlying, GridFittedBottom(island))

    for grid in (periodic, bounded, immersed)
        @info "  Comparing skipping and non-skipping solvers on $(summary(grid))..."
        full, skip, same = compare_skipping(grid)

        # Sanity check: the ice moved and part of the domain is ice free
        @test maximum(abs, interior(full.velocities.u)) > 0
        @test any(iszero, interior(full.ice_thickness))

        for (name, s) in same
            @test s
        end

        # Once the state has settled, the first substep only computes the iced points
        maps = skip.dynamics.solver.ice_free_cells.velocity
        @test maps.host_counts[1] < length(maps.first_substep)
    end

    @info "  Testing velocities set over open water..."
    _, _, same = compare_skipping(periodic; perturb_velocities = true)
    @test all(values(same))

    @info "  Testing the limits of no ice and full ice cover..."
    for ice in (:none, :full)
        _, _, same = compare_skipping(periodic; ice, Nt = 3)
        @test all(values(same))
    end

    @info "  Testing a rheology without stress kernels..."
    _, _, same = compare_skipping(periodic; rheology = ViscousRheology(ν = 1000), Nt = 3)
    @test all(values(same))
end
