using Test
using Oceananigans
using Oceananigans.Advection: EnergyConserving, EnstrophyConserving
using Oceananigans.BoundaryConditions: FluxBoundaryCondition, fill_halo_regions!
using Oceananigans.Coriolis: ActiveWeightedEnergyConserving, ActiveWeightedEnstrophyConserving,
                             TriadScheme, x_f_cross_U, y_f_cross_U
using Oceananigans.Fields: ConstantField
using Oceananigans.Grids: MutableVerticalDiscretization
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryCondition
using ClimaSeaIce
using ClimaSeaIce.Rheologies: immersed_∂ⱼ_σ₁ⱼ, immersed_∂ⱼ_σ₂ⱼ
using ClimaSeaIce.SeaIceDynamics: x_f_cross_U_2D, y_f_cross_U_2D

const Nx, Ny, Nz = 12, 12, 4
const zstretched = [-100, -60, -30, -10, 0]

coriolis_schemes = (EnstrophyConserving(),
                    EnergyConserving(),
                    ActiveWeightedEnstrophyConserving(),
                    ActiveWeightedEnergyConserving(),
                    TriadScheme())

underlying_grid(z) = LatitudeLongitudeGrid(size = (Nx, Ny, Nz), longitude = (0, 30), latitude = (60, 82),
                                           z = z, halo = (4, 4, 4), topology = (Bounded, Bounded, Bounded))

# An island that breaches the surface layer, so that the top layer carries immersed cells.
bathymetry(λ, φ) = -100 + 130 * exp(-((λ - 15)^2 + (φ - 70)^2) / 20)

immersed_grid(z) = ImmersedBoundaryGrid(underlying_grid(z), GridFittedBottom(bathymetry))

# A free surface of ±1 m on a 100 m column, written into the vertical scaling of a mutable grid.
function set_surface_scaling!(grid; H = 100, ηmax = 1)
    for σ in (grid.z.σᶜᶜⁿ, grid.z.σᶠᶜⁿ, grid.z.σᶜᶠⁿ, grid.z.σᶠᶠⁿ, grid.z.σᶜᶜ⁻)
        parent(σ) .= 1
    end

    η(i, j) = ηmax * sin(2π * i / Nx) * cos(2π * j / Ny)

    for i in axes(parent(grid.z.σᶜᶜⁿ), 1), j in axes(parent(grid.z.σᶜᶜⁿ), 2)
        parent(grid.z.σᶜᶜⁿ)[i, j, 1] = 1 + η(i, j) / H
        parent(grid.z.σᶠᶜⁿ)[i, j, 1] = 1 + η(i - 1/2, j) / H
        parent(grid.z.σᶜᶠⁿ)[i, j, 1] = 1 + η(i, j - 1/2) / H
        parent(grid.z.σᶠᶠⁿ)[i, j, 1] = 1 + η(i - 1/2, j - 1/2) / H
        parent(grid.z.σᶜᶜ⁻)[i, j, 1] = parent(grid.z.σᶜᶜⁿ)[i, j, 1]
    end

    parent(grid.z.ηⁿ) .= 0

    return grid
end

relative_difference(a, b) = maximum(abs, a .- b) / maximum(abs, a)

function sea_ice_velocities(grid)
    u = Field{Face, Center, Nothing}(grid)
    v = Field{Center, Face, Nothing}(grid)
    set!(u, (λ, φ) -> 0.1 * sind(12λ) * cosd(9φ))
    set!(v, (λ, φ) -> 0.07 * cosd(7λ) * sind(11φ))
    fill_halo_regions!(u)
    fill_halo_regions!(v)
    return (; u, v)
end

@testset "Two-dimensional Coriolis" begin
    @info "Testing that the sea ice Coriolis acceleration carries no vertical spacing..."

    # On a grid whose Δz is horizontally uniform, the two-dimensional operators must reproduce the
    # area-weighted ones of Oceananigans exactly.
    for grid in (underlying_grid(zstretched), immersed_grid(zstretched))
        U = sea_ice_velocities(grid)
        kᴺ = size(grid, 3)

        for scheme in coriolis_schemes
            coriolis = HydrostaticSphericalCoriolis(; scheme)

            for i in 2:Nx-1, j in 2:Ny-1
                @test x_f_cross_U_2D(i, j, kᴺ, grid, coriolis, U) ≈ x_f_cross_U(i, j, kᴺ, grid, coriolis, U)
                @test y_f_cross_U_2D(i, j, kᴺ, grid, coriolis, U) ≈ y_f_cross_U(i, j, kᴺ, grid, coriolis, U)
            end
        end
    end

    # A moving grid must leave the acceleration untouched.
    for build_grid in (underlying_grid, immersed_grid)
        static = build_grid(zstretched)
        moving = set_surface_scaling!(build_grid(MutableVerticalDiscretization(zstretched)))

        Us = sea_ice_velocities(static)
        Um = sea_ice_velocities(moving)
        kᴺ = size(static, 3)

        for scheme in coriolis_schemes
            coriolis = HydrostaticSphericalCoriolis(; scheme)

            for i in 2:Nx-1, j in 2:Ny-1
                @test x_f_cross_U_2D(i, j, kᴺ, moving, coriolis, Um) == x_f_cross_U_2D(i, j, kᴺ, static, coriolis, Us)
                @test y_f_cross_U_2D(i, j, kᴺ, moving, coriolis, Um) == y_f_cross_U_2D(i, j, kᴺ, static, coriolis, Us)
            end
        end
    end
end

@testset "Immersed stress divergence" begin
    @info "Testing that the immersed stress divergence carries no vertical spacing..."

    wall_stress = ImmersedBoundaryCondition(west  = FluxBoundaryCondition(0.3),
                                            east  = FluxBoundaryCondition(0.3),
                                            south = FluxBoundaryCondition(0.2),
                                            north = FluxBoundaryCondition(0.2))

    static = immersed_grid(zstretched)
    moving = set_surface_scaling!(immersed_grid(MutableVerticalDiscretization(zstretched)))

    kᴺ = size(static, 3)
    rheology = ViscousRheology(ν = 0)
    clock = nothing
    fields = sea_ice_velocities(static)

    σ₁ = [immersed_∂ⱼ_σ₁ⱼ(i, j, kᴺ, static, wall_stress, rheology, clock, fields) for i in 2:Nx-1, j in 2:Ny-1]
    σ₂ = [immersed_∂ⱼ_σ₂ⱼ(i, j, kᴺ, static, wall_stress, rheology, clock, fields) for i in 2:Nx-1, j in 2:Ny-1]

    # The comparison below is only worth making where the wall actually exerts a stress.
    @test any(!iszero, σ₁)
    @test any(!iszero, σ₂)

    for i in 2:Nx-1, j in 2:Ny-1
        @test immersed_∂ⱼ_σ₁ⱼ(i, j, kᴺ, moving, wall_stress, rheology, clock, fields) ==
              immersed_∂ⱼ_σ₁ⱼ(i, j, kᴺ, static, wall_stress, rheology, clock, fields)
        @test immersed_∂ⱼ_σ₂ⱼ(i, j, kᴺ, moving, wall_stress, rheology, clock, fields) ==
              immersed_∂ⱼ_σ₂ⱼ(i, j, kᴺ, static, wall_stress, rheology, clock, fields)
    end
end

@testset "Time stepping on a moving vertical grid" begin
    @info "Testing that a sea ice time step is insensitive to the vertical grid..."

    function run_sea_ice_model(grid; timestepper, advection, nsteps = 5, Δt = 300)
        τₒ = SemiImplicitStress(uₑ = ConstantField(0.15), vₑ = ConstantField(-0.05))

        dynamics = SeaIceMomentumEquation(grid;
                                          coriolis = HydrostaticSphericalCoriolis(),
                                          bottom_momentum_stress = τₒ,
                                          rheology = ElastoViscoPlasticRheology(),
                                          solver = SplitExplicitSolver(grid; substeps = 20))

        model = SeaIceModel(grid; dynamics, advection, timestepper)

        set!(model, h = (λ, φ) -> 1 + 0.5 * sind(18λ) * cosd(24φ),
                    ℵ = (λ, φ) -> 0.6 + 0.3 * sind(24φ),
                    u = 0.01, v = -0.02)

        for _ in 1:nsteps
            time_step!(model, Δt)
        end

        return map(Array, (h = interior(model.ice_thickness, :, :, 1),
                           ℵ = interior(model.ice_concentration, :, :, 1),
                           u = interior(model.velocities.u, :, :, 1),
                           v = interior(model.velocities.v, :, :, 1)))
    end

    for build_grid in (underlying_grid, immersed_grid),
        (timestepper, advection) in ((:SplitRungeKutta3, WENO(order=5)),
                                     (:ForwardEuler, IncrementalRemapping()))

        static = build_grid(zstretched)
        moving = set_surface_scaling!(build_grid(MutableVerticalDiscretization(zstretched)))

        Ψs = run_sea_ice_model(static; timestepper, advection)
        Ψm = run_sea_ice_model(moving; timestepper, advection)

        for name in keys(Ψs)
            @test relative_difference(Ψs[name], Ψm[name]) < 1e-12
        end
    end
end
