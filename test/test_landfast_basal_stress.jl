using Test
using Oceananigans
using ClimaSeaIce

using ClimaSeaIce.SeaIceDynamics: basal_τx_coefficient, basal_τy_coefficient,
                                  materialize_basal_stress, basal_stress_magnitude
using Oceananigans.Fields: interior
using Oceananigans.Grids: Face, Center

# A shelf shallow enough to ground keels next to water too deep to ground anything. The vertical grid
# must resolve the shelf: with cells thicker than the shelf depth the whole top cell is immersed and
# the column depth collapses to zero.
function shelf_grid(arch = CPU(); Nx = 8, Ny = 8, Nz = 8)
    underlying = RectilinearGrid(arch; size = (Nx, Ny, Nz), x = (0, 8e4), y = (0, 8e4), z = (-40, 0),
                                 topology = (Bounded, Bounded, Bounded))
    bottom = [i <= Nx ÷ 2 ? -20.0 : -40.0 for i in 1:Nx, j in 1:Ny]
    return ImmersedBoundaryGrid(underlying, GridFittedBottom(bottom))
end

ice_fields(grid; h = 3.0, ℵ = 1.0, u = 0.1, v = 0.0) = begin
    hf = Field{Center, Center, Nothing}(grid); set!(hf, h)
    ℵf = Field{Center, Center, Nothing}(grid); set!(ℵf, ℵ)
    ρf = Field{Center, Center, Nothing}(grid); set!(ρf, 900.0)
    uf = Field{Face,   Center, Nothing}(grid); set!(uf, u)
    vf = Field{Center, Face,   Nothing}(grid); set!(vf, v)
    (; h = hf, ℵ = ℵf, ρ = ρf, u = uf, v = vf)
end

@testset "LandfastBasalStress construction" begin
    b = LandfastBasalStress()
    @test b.water_depth isa Nothing
    @test b.critical_thickness_parameter == 8
    @test b.stress_parameter == 15
    @test b.concentration_hardening == 20
    @test b.minimum_speed == 5e-5
    @test b.maximum_water_depth == 30

    grid = shelf_grid()
    bm = materialize_basal_stress(b, grid)
    H = interior(bm.water_depth)[:, :, 1]
    @test all(H[1:4, :] .≈ 20)
    @test all(H[5:8, :] .≈ 40)

    @test materialize_basal_stress(nothing, grid) isa Nothing
end

@testset "Grounding criterion" begin
    grid = shelf_grid()
    bm = materialize_basal_stress(LandfastBasalStress(), grid)
    kᴺ = size(grid, 3)

    fields = ice_fields(grid)   # h = 3 m, ℵ = 1, u = 0.1 m/s
    shelf = basal_τx_coefficient(2, 4, kᴺ, grid, bm, fields)
    deep  = basal_τx_coefficient(7, 4, kᴺ, grid, bm, fields)

    # hᶜ = H ℵ / k₁ = 20/8 = 2.5 m, so 3 m of ice grounds on the shelf.
    @test shelf ≈ 15 * (3.0 - 20/8) / (0.1 + 5e-5)
    @test deep == 0                       # 40 m exceeds maximum_water_depth

    # τᵇ saturates at k₂ δh: it can arrest ice but never accelerate it.
    @test shelf * 0.1 < 15 * (3.0 - 20/8)

    # Thin ice cannot ground even on the shelf.
    @test basal_τx_coefficient(2, 4, kᴺ, grid, materialize_basal_stress(LandfastBasalStress(), grid),
                               ice_fields(grid; h = 1.0)) == 0

    # Unconsolidated ice is released by the concentration hardening.
    loose = basal_τx_coefficient(2, 4, kᴺ, grid, bm, ice_fields(grid; ℵ = 0.5))
    @test 0 < loose < shelf

    # Slower ice feels a larger coefficient, so the stress saturates rather than the drag vanishing.
    @test basal_τx_coefficient(2, 4, kᴺ, grid, bm, ice_fields(grid; u = 0.01)) > shelf

    @test basal_τx_coefficient(2, 4, kᴺ, grid, nothing, fields) == 0
    @test basal_τy_coefficient(2, 4, kᴺ, grid, nothing, fields) == 0
end

@testset "Grounding is a cell-centred property" begin
    grid = shelf_grid()
    bm = materialize_basal_stress(LandfastBasalStress(), grid)
    kᴺ = size(grid, 3)
    fields = ice_fields(grid)

    # The magnitude is formed at centres, so the shelf cell carries the full grounding, and the
    # velocity point straddling the shelf break gets the average of its two neighbours.
    mˢ = basal_stress_magnitude(4, 4, kᴺ, grid, bm, fields)
    mᵈ = basal_stress_magnitude(5, 4, kᴺ, grid, bm, fields)
    @test mˢ ≈ 15 * (3.0 - 20/8)
    @test mᵈ == 0
    @test basal_τx_coefficient(5, 4, kᴺ, grid, bm, fields) ≈ (mˢ + mᵈ) / 2 / (0.1 + 5e-5)
end

@testset "SeaIceMomentumEquation carries the basal stress" begin
    grid = shelf_grid()

    plain = SeaIceMomentumEquation(grid)
    @test plain.basal_stress isa Nothing

    dyn = SeaIceMomentumEquation(grid; basal_stress = LandfastBasalStress())
    @test dyn.basal_stress isa LandfastBasalStress
    @test dyn.basal_stress.water_depth isa Field       # materialized on the grid

    # It is held apart from the stresses the water column exerts, so a coupler reading those never
    # sees the part of the drag that the sea floor carries.
    @test !(:basal in propertynames(dyn.external_momentum_stresses))
end

function drift_after_a_step(; grounded, Δt = 600, nsteps = 20)
    grid = shelf_grid()
    basal_stress = grounded ? LandfastBasalStress() : nothing
    dynamics = SeaIceMomentumEquation(grid; coriolis = nothing, basal_stress,
                                      rheology = ElastoViscoPlasticRheology(),
                                      solver = SplitExplicitSolver(grid; substeps = 50))
    model = SeaIceModel(grid; dynamics, advection = nothing)
    set!(model, h = 3.0, ℵ = 1.0)
    set!(model.velocities.u, 0.1)
    for _ in 1:nsteps
        time_step!(model, Δt)
    end
    return Array(interior(model.velocities.u))[:, :, 1]
end

@testset "Grounded ice is arrested, deep ice is not" begin
    uᶠ = drift_after_a_step(grounded = false)
    uᵍ = drift_after_a_step(grounded = true)

    @test all(isfinite, uᶠ)
    @test all(isfinite, uᵍ)

    # Columns 1:4 sit on the 20 m shelf and ground; 6:8 are in 40 m water and cannot. The deep ice
    # is not perfectly untouched: ice is a continuum, so the arrested belt drags on its neighbours
    # through the internal stress. What must hold is the contrast, not decoupling.
    shelf_slowdown = 1 - maximum(abs, uᵍ[2:4, :]) / maximum(abs, uᶠ[2:4, :])
    deep_slowdown  = 1 - abs(uᵍ[7, 4]) / abs(uᶠ[7, 4])

    @test shelf_slowdown > 0.5
    @test deep_slowdown < 0.05
    @test shelf_slowdown > 10 * deep_slowdown
end
