using Test
using Oceananigans
using ClimaSeaIce

using ClimaSeaIce.Rheologies: u_lateral_boundary_condition, v_lateral_boundary_condition,
                              lateral_boundary_condition, wall_shear_stress, strain_rate_xy,
                              _ice_stress_uy, _ice_stress_vx, _ice_stress_ux, _ice_stress_vy
using Oceananigans.ImmersedBoundaries: immersed_peripheral_node, immersed_inactive_node
using Oceananigans.Operators: Δyᶠᶠᶜ
using Oceananigans.Grids: Face, Center, peripheral_node

# A channel walled in y, so every interior column has an immersed face to its north or south.
function channel_grid(arch = CPU(); Nx = 8, Ny = 8, Nz = 1)
    underlying = RectilinearGrid(arch; size = (Nx, Ny, Nz), x = (0, 8e4), y = (0, 8e4), z = (-10, 0),
                                 topology = (Periodic, Bounded, Bounded))
    bottom = [(j == 1 || j == Ny) ? 0.0 : -10.0 for i in 1:Nx, j in 1:Ny]
    return ImmersedBoundaryGrid(underlying, GridFittedBottom(bottom))
end

# No slip is requested exactly as on a domain boundary: a zero-valued boundary condition.
function velocity_bcs(grid; no_slip)
    immersed = no_slip ? ValueBoundaryCondition(0) : nothing
    u = FieldBoundaryConditions(grid, (Face(), Center(), nothing); immersed)
    v = FieldBoundaryConditions(grid, (Center(), Face(), nothing); immersed)
    return (; u, v)
end

@testset "A zero-valued immersed BC means no slip" begin
    grid = channel_grid()

    free = velocity_bcs(grid; no_slip = false)
    nsl  = velocity_bcs(grid; no_slip = true)

    # Oceananigans puts the Value on the sides where each component is tangential: south/north
    # for u, west/east for v.
    @test u_lateral_boundary_condition(free.u.immersed) isa FreeSlip
    @test v_lateral_boundary_condition(free.v.immersed) isa FreeSlip
    @test u_lateral_boundary_condition(nsl.u.immersed) isa NoSlip
    @test v_lateral_boundary_condition(nsl.v.immersed) isa NoSlip
    @test u_lateral_boundary_condition(nothing) isa FreeSlip
    @test v_lateral_boundary_condition(nothing) isa FreeSlip

    # The shear strain rate takes either component asking for no slip.
    @test lateral_boundary_condition(nsl.u.immersed, nsl.v.immersed) isa NoSlip
    @test lateral_boundary_condition(nsl.u.immersed, free.v.immersed) isa NoSlip
    @test lateral_boundary_condition(free.u.immersed, free.v.immersed) isa FreeSlip

    # A flux boundary condition is a user stress, not a slip condition.
    flux = FieldBoundaryConditions(grid, (Face(), Center(), nothing);
                                   immersed = FluxBoundaryCondition(0))
    @test u_lateral_boundary_condition(flux.immersed) isa FreeSlip

    # The rheology no longer carries the condition.
    @test !hasproperty(ElastoViscoPlasticRheology(), :lateral_boundary_condition)
end

# Differencing against the zero in a land cell puts the wall half a cell inside the land. Imposing
# u = 0 on the face itself needs the reflected value -u, which doubles the gradient.
@testset "No-slip doubles the wall shear strain rate" begin
    grid = channel_grid()
    kᴺ = size(grid, 3)

    u = Field{Face,   Center, Nothing}(grid)
    v = Field{Center, Face,   Nothing}(grid)
    set!(u, 0.1)
    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        if peripheral_node(i, j, kᴺ, grid, Face(), Center(), Center())
            u[i, j, 1] = 0
        end
    end

    i, jʷ, jᵢ = 4, 2, 4
    Δy = Δyᶠᶠᶜ(i, jʷ, kᴺ, grid)
    uᶜ = u[i, jʷ, 1]

    @test u[i, jʷ - 1, 1] == 0                       # the land cell is masked, not reflected

    ϵ̇ᶠ = strain_rate_xy(i, jʷ, kᴺ, grid, u, v, FreeSlip())
    ϵ̇ⁿ = strain_rate_xy(i, jʷ, kᴺ, grid, u, v, NoSlip())

    @test ϵ̇ᶠ ≈ uᶜ / (2Δy)                            # wall effectively half a cell inside the land
    @test ϵ̇ⁿ ≈ uᶜ / Δy                               # u = 0 on the face, as no slip requires
    @test ϵ̇ⁿ ≈ 2ϵ̇ᶠ

    @test strain_rate_xy(i, jᵢ, kᴺ, grid, u, v, NoSlip()) ==
          strain_rate_xy(i, jᵢ, kᴺ, grid, u, v, FreeSlip())
end

# A wall normal to x: v is tangential to it and must be reflected, u is the normal component whose
# zero is exact impermeability and must NOT be.
function strip_grid(; Nx = 8, Ny = 8, Nz = 1)
    underlying = RectilinearGrid(size = (Nx, Ny, Nz), x = (0, 8e4), y = (0, 8e4), z = (-10, 0),
                                 topology = (Periodic, Periodic, Bounded))
    bottom = [(i == 4 || i == 5) ? 0.0 : -10.0 for i in 1:Nx, j in 1:Ny]
    return ImmersedBoundaryGrid(underlying, GridFittedBottom(bottom))
end

@testset "Impermeable faces are not mistaken for land" begin
    grid = strip_grid()
    kᴺ = size(grid, 3)
    i, j = 4, 4      # the western face of the solid strip

    # The u-node there is peripheral — u is zero by impermeability — but it is not inside the solid.
    @test immersed_peripheral_node(i, j, kᴺ, grid, Face(), Center(), Center())
    @test !immersed_inactive_node(i, j, kᴺ, grid, Face(), Center(), Center())

    # The v-node one column in *is* buried in the solid, so its zero stands in for a reflection.
    @test immersed_inactive_node(i, j, kᴺ, grid, Center(), Face(), Center())

    u = Field{Face,   Center, Nothing}(grid)
    v = Field{Center, Face,   Nothing}(grid)
    set!(v, 0.1)
    for jj in 1:size(grid, 2), ii in 1:size(grid, 1)
        peripheral_node(ii, jj, kᴺ, grid, Center(), Face(), Center()) && (v[ii, jj, 1] = 0)
    end

    Δx = Oceananigans.Operators.Δxᶠᶠᶜ(i, j, kᴺ, grid)
    vᶜ = v[i - 1, j, 1]

    ϵ̇ᶠ = strain_rate_xy(i, j, kᴺ, grid, u, v, FreeSlip())
    ϵ̇ⁿ = strain_rate_xy(i, j, kᴺ, grid, u, v, NoSlip())

    @test ϵ̇ⁿ ≈ 2ϵ̇ᶠ            # the tangential v gradient is reflected
    @test ϵ̇ᶠ ≈ -vᶜ / (2Δx)
    @test ϵ̇ⁿ ≈ -vᶜ / Δx
end

@testset "Only the shear stress feels the wall" begin
    grid = channel_grid()
    kᴺ = size(grid, 3)
    i, jʷ, jᵢ = 4, 2, 4

    @test immersed_peripheral_node(i, jʷ, kᴺ, grid, Face(), Face(), Center())
    @test !immersed_peripheral_node(i, jᵢ, kᴺ, grid, Face(), Face(), Center())

    σ = 12345.0
    @test wall_shear_stress(FreeSlip(), i, jʷ, kᴺ, grid, σ) == 0
    @test wall_shear_stress(NoSlip(),   i, jʷ, kᴺ, grid, σ) == σ
    @test wall_shear_stress(FreeSlip(), i, jᵢ, kᴺ, grid, σ) == σ
    @test wall_shear_stress(NoSlip(),   i, jᵢ, kᴺ, grid, σ) == σ

    σ₁₁ = Field{Center, Center, Nothing}(grid); set!(σ₁₁, 3.0)
    σ₂₂ = Field{Center, Center, Nothing}(grid); set!(σ₂₂, 3.0)
    σ₁₂ = Field{Face,   Face,   Nothing}(grid); set!(σ₁₂, 7.0)
    fields = (; σ₁₁, σ₂₂, σ₁₂)
    clock = Clock(grid)
    r = ElastoViscoPlasticRheology()

    @test _ice_stress_uy(i, jʷ, kᴺ, grid, r, clock, fields, FreeSlip()) == 0
    @test _ice_stress_uy(i, jʷ, kᴺ, grid, r, clock, fields, NoSlip())   == 7
    @test _ice_stress_vx(i, jʷ, kᴺ, grid, r, clock, fields, FreeSlip()) == 0
    @test _ice_stress_vx(i, jʷ, kᴺ, grid, r, clock, fields, NoSlip())   == 7

    # Normal stresses are masked identically under both conditions.
    for lbc in (FreeSlip(), NoSlip())
        @test _ice_stress_ux(i, 1, kᴺ, grid, r, clock, fields) == 0     # dry row
        @test _ice_stress_vy(i, 1, kᴺ, grid, r, clock, fields) == 0
        @test _ice_stress_ux(i, jᵢ, kᴺ, grid, r, clock, fields) == 3
        @test _ice_stress_vy(i, jᵢ, kᴺ, grid, r, clock, fields) == 3
    end
end

function drift_after_a_step(; no_slip, substeps = 50, Δt = 100)
    grid = channel_grid()
    dynamics = SeaIceMomentumEquation(grid; coriolis = nothing,
                                      rheology = ElastoViscoPlasticRheology(),
                                      solver = SplitExplicitSolver(grid; substeps))
    model = SeaIceModel(grid; dynamics, advection = nothing,
                        boundary_conditions = velocity_bcs(grid; no_slip))

    set!(model, h = 1.0, ℵ = 1.0)
    set!(model.velocities.u, 0.1)
    time_step!(model, Δt)

    return Array(interior(model.velocities.u))[:, :, 1]
end

@testset "No-slip decelerates ice along a coast" begin
    uᶠ = drift_after_a_step(no_slip = false)
    uⁿ = drift_after_a_step(no_slip = true)

    @test all(isfinite, uᶠ)
    @test all(isfinite, uⁿ)

    # The wall can only remove momentum, and it must remove some in the row next to the coast.
    @test all(abs.(uⁿ) .<= abs.(uᶠ) .+ 1e-12)
    @test maximum(abs.(uᶠ[:, 2]) .- abs.(uⁿ[:, 2])) > 0
end
