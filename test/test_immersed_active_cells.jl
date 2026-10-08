using Test
using Oceananigans
using Oceananigans.Units: hour
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: φnode
using Oceananigans.OrthogonalSphericalShellGrids: TripolarGrid
using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics: SplitExplicitSolver, SemiImplicitStress, SeaIceMomentumEquation
using ClimaSeaIce.Rheologies: ElastoViscoPlasticRheology

# Sea ice can only exist poleward of `latitude_cutoff`: everything equatorward of it is immersed,
# as is some land, so that the split-explicit solver skips most of the domain.
const latitude_cutoff = 40

function band_bottom(λ, φ)
    land = (abs(φ) < latitude_cutoff) | ((60 < φ < 70) & (100 < λ < 140))
    return land ? 0.0 : -10.0
end

latitude_longitude_grid(; active_cells_map) =
    ImmersedBoundaryGrid(LatitudeLongitudeGrid(size = (72, 40, 1), longitude = (0, 360), latitude = (-80, 80),
                                               z = (-10, 0), halo = (7, 7, 7)),
                         GridFittedBottom(band_bottom); active_cells_map)

function tripolar_grid(; active_cells_map)
    underlying = TripolarGrid(CPU(); size = (60, 60, 1), southernmost_latitude = -75, halo = (7, 7, 7), z = (-10, 0))
    λp = underlying.conformal_mapping.first_pole_longitude
    φp = underlying.conformal_mapping.north_poles_latitude

    # Land around the two northern poles, as on an ORCA grid
    pole(λ, φ) = (abs(φp - φ) < 5) & ((abs(λ - λp) < 5) | (abs(λ - λp - 180) < 5) | (abs(λ - λp - 360) < 5))
    bottom(λ, φ) = pole(λ, φ) ? 0.0 : band_bottom(λ, φ)

    return ImmersedBoundaryGrid(underlying, GridFittedBottom(bottom); active_cells_map)
end

function immersed_band_model(grid; rheology = ElastoViscoPlasticRheology(), no_slip = false)
    uₒ = XFaceField(grid)
    vₒ = YFaceField(grid)
    set!(uₒ, (λ, φ, z...) ->  0.2 * cosd(3λ))
    set!(vₒ, (λ, φ, z...) -> -0.2 * sind(2λ))
    fill_halo_regions!((uₒ, vₒ))

    dynamics = SeaIceMomentumEquation(grid; coriolis = HydrostaticSphericalCoriolis(),
                                      bottom_momentum_stress = SemiImplicitStress(uₑ = uₒ, vₑ = vₒ),
                                      rheology, solver = SplitExplicitSolver(grid; substeps = 50))

    immersed = no_slip ? ValueBoundaryCondition(0) : nothing
    boundary_conditions = (u = FieldBoundaryConditions(grid, (Face(), Center(), nothing); immersed),
                           v = FieldBoundaryConditions(grid, (Center(), Face(), nothing); immersed))

    model = SeaIceModel(grid; dynamics, boundary_conditions,
                        advection = ClimaSeaIce.IncrementalRemapping(), timestepper = :ForwardEuler)

    ice(φ) = abs(φ) > latitude_cutoff + 10
    set!(model, h = (λ, φ, z...) -> ice(φ) * (1 + 0.5 * sind(2λ)),
                ℵ = (λ, φ, z...) -> ice(φ) * 0.9)
    return model
end

# Run the same model on the same immersed grid, with and without the active-columns map
function compare_active_cells(build_grid; Nt = 6, kw...)
    full = immersed_band_model(build_grid(active_cells_map = false); kw...)
    skip = immersed_band_model(build_grid(active_cells_map = true);  kw...)

    for _ in 1:Nt
        time_step!(full, 1hour)
        time_step!(skip, 1hour)
    end

    names = (:u, :v, :h, :ℵ)
    # `isequal` so that NaNs in the same places of both (e.g. unused tripolar halo corners) count as equal
    same = Dict(name => isequal(parent(fields(full)[name]), parent(fields(skip)[name])) for name in names)

    aux_full = full.dynamics.auxiliaries.fields
    aux_skip = skip.dynamics.auxiliaries.fields
    for name in keys(aux_full)
        name == :Δ && continue # diagnostic, never computed in immersed columns

        # The viscosities are only read where they were just computed, so they are not computed at the
        # skipped points either (where the full kernels can leave NaNs): compare them at the computed points
        same[name] = if name in (:ζᶜᶜᶜ, :ζᶠᶠᶜ)
            computed = Array(skip.dynamics.solver.active_cells.lists.stress)
            all(isequal(aux_full[name][i, j, 1], aux_skip[name][i, j, 1]) for (i, j) in computed)
        else
            isequal(parent(aux_full[name]), parent(aux_skip[name]))
        end
    end

    return full, skip, same
end

@testset "Split-explicit solver on an immersed grid with an active-columns map" begin
    for (name, build_grid) in (("LatitudeLongitudeGrid", latitude_longitude_grid), ("TripolarGrid", tripolar_grid))
        @info "  Comparing the solver with and without the active-columns map on a $name..."
        full, skip, same = compare_active_cells(build_grid)

        # Only the grid with the map skips immersed columns, and it skips most of the domain
        @test isnothing(full.dynamics.solver.active_cells)
        lists = skip.dynamics.solver.active_cells.lists
        @test length(lists.velocity) < 0.6 * size(skip.grid, 1) * size(skip.grid, 2)
        @test !isnothing(lists.stress)

        # Sanity checks: the ice moved, and there is none in the immersed band
        @test maximum(abs, interior(full.velocities.u)) > 0
        Nx, Ny, _ = size(skip.grid)
        h = interior(skip.ice_thickness, :, :, 1)
        @test all(iszero(h[i, j]) for i in 1:Nx, j in 1:Ny if abs(φnode(i, j, 1, skip.grid, Center(), Center(), Center())) < latitude_cutoff)

        for (field, s) in same
            s || @warn "  $field differs on the $name"
            @test s
        end
    end

    @info "  Comparing with no-slip walls..."
    for build_grid in (latitude_longitude_grid, tripolar_grid)
        _, _, same = compare_active_cells(build_grid; no_slip = true, Nt = 3)
        @test all(values(same))
    end

    @info "  Comparing a rheology without stress kernels..."
    _, skip, same = compare_active_cells(latitude_longitude_grid; rheology = ViscousRheology(ν = 1000), Nt = 3)
    @test isnothing(skip.dynamics.solver.active_cells.lists.stress)
    @test all(values(same))
end
