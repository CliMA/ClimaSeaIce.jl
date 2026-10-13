using Test
using Oceananigans
using Oceananigans.Units: hour
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: φnode, surface_kernel_parameters
using Oceananigans.OrthogonalSphericalShellGrids: TripolarGrid
using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics: SplitExplicitSolver, SemiImplicitStress, SeaIceMomentumEquation
using ClimaSeaIce.Rheologies: ElastoViscoPlasticRheology

# Sea ice can only exist poleward of `latitude_cutoff`: everything equatorward of it is immersed, as is
# some land (and, on the tripolar grid, land around the northern poles, as on an ORCA grid)
const latitude_cutoff = 40

land(λ, φ) = (abs(φ) < latitude_cutoff) | ((60 < φ < 70) & (100 < λ < 140))

latitude_longitude_bottom(underlying) = (λ, φ) -> land(λ, φ) ? 0.0 : -10.0

function tripolar_bottom(underlying)
    λp = underlying.conformal_mapping.first_pole_longitude
    φp = underlying.conformal_mapping.north_poles_latitude
    pole(λ, φ) = (abs(φp - φ) < 5) & any(abs(λ - λp - Δλ) < 5 for Δλ in (0, 180, 360))
    return (λ, φ) -> land(λ, φ) | pole(λ, φ) ? 0.0 : -10.0
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
    set!(model, h = (λ, φ, z...) -> ice(φ) * (1 + 0.5 * sind(2λ)), ℵ = (λ, φ, z...) -> ice(φ) * 0.9)
    return model
end

latitude_longitude = LatitudeLongitudeGrid(size = (72, 40, 1), longitude = (0, 360), latitude = (-80, 80),
                                           z = (-10, 0), halo = (7, 7, 7))
tripolar = TripolarGrid(CPU(); size = (60, 60, 1), southernmost_latitude = -75, halo = (7, 7, 7), z = (-10, 0))

cases = (("LatitudeLongitudeGrid",           latitude_longitude, latitude_longitude_bottom, (;)),
         ("TripolarGrid",                    tripolar,           tripolar_bottom,           (;)),
         ("LatitudeLongitudeGrid, no slip",  latitude_longitude, latitude_longitude_bottom, (; no_slip = true)),
         ("TripolarGrid, no slip",           tripolar,           tripolar_bottom,           (; no_slip = true)),
         ("LatitudeLongitudeGrid, viscous",  latitude_longitude, latitude_longitude_bottom, (; rheology = ViscousRheology(ν = 1000))))

# The same model on the same immersed grid, with and without the active-columns map, must agree bitwise
@testset "Split-explicit solver on an immersed grid with an active-columns map" begin
    for (name, underlying, bottom, kw) in cases
        @testset "$name" begin
            bottom_height = bottom(underlying)
            full = immersed_band_model(ImmersedBoundaryGrid(underlying, GridFittedBottom(bottom_height); active_cells_map = false); kw...)
            skip = immersed_band_model(ImmersedBoundaryGrid(underlying, GridFittedBottom(bottom_height); active_cells_map = true);  kw...)

            for _ in 1:4
                time_step!(full, 1hour)
                time_step!(skip, 1hour)
            end

            # The map skips most of the domain, the ice moved, and there is none in the immersed band
            Nx, Ny, _ = size(skip.grid)
            @test length(skip.grid.active_z_columns) < 0.6 * Nx * Ny
            @test maximum(abs, interior(full.velocities.u)) > 0
            @test all(iszero(skip.ice_thickness[i, j, 1]) for i in 1:Nx, j in 1:Ny
                      if abs(φnode(i, j, 1, skip.grid, Center(), Center(), Center())) < latitude_cutoff)

            # `isequal` so that NaNs in the same places of both (e.g. unused tripolar halo corners) count as equal
            for field in (:u, :v, :h, :ℵ)
                @test isequal(parent(fields(full)[field]), parent(fields(skip)[field]))
            end

            # Δ is diagnostic and never computed in immersed columns. The viscosities are only read where they
            # were just computed, so they are not computed at the skipped points either (where the full kernels
            # can leave NaNs): compare them at the computed points
            aux_full = full.dynamics.auxiliaries.fields
            aux_skip = skip.dynamics.auxiliaries.fields
            computed = Array(surface_kernel_parameters(skip.grid))
            for field in keys(aux_full)
                if field in (:ζᶜᶜᶜ, :ζᶠᶠᶜ)
                    @test all(isequal(aux_full[field][i, j, 1], aux_skip[field][i, j, 1]) for (i, j) in computed)
                elseif field != :Δ
                    @test isequal(parent(aux_full[field]), parent(aux_skip[field]))
                end
            end
        end
    end
end
