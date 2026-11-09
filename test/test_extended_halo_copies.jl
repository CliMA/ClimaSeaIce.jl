using Test
using Oceananigans
using ClimaSeaIce

using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: halo_size
using Oceananigans.ImmersedBoundaries: has_active_cells_map
using ClimaSeaIce.SeaIceDynamics: extended_external_variable, extended_sea_ice_grid, refresh_and_fill_external_field!

# On a distributed grid, the split-explicit (mEVP) solver extends the halos of the sea-ice velocity grid to
# 2 substeps + 3 cells and copies the external (ocean) velocities onto it. The ocean surface velocity is a k = Nz view
# of a three-dimensional field: its extended copy must stay a two-dimensional slab, and the extended immersed grid
# must not carry a three-dimensional active-cells map, since sea-ice fields are two-dimensional.
@testset "EVP-extended copies stay two-dimensional" begin
    @info "Testing the memory of the EVP-extended copies of external fields..."

    Nx, Ny, Nz = 32, 32, 10
    underlying_grid = RectilinearGrid(size = (Nx, Ny, Nz), halo = (4, 4, 3), x = (0, 1e5), y = (0, 1e5), z = (-1000, 0),
                                      topology = (Periodic, Bounded, Bounded))
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom((x, y) -> y > 8e4 ? 0 : -500); active_cells_map = true)

    halo = (27, 27, 3) # 2 × 12 substeps + 3
    extended_grid = extended_sea_ice_grid(halo, grid)
    @test halo_size(extended_grid) == halo
    @test extended_grid isa ImmersedBoundaryGrid
    @test has_active_cells_map(grid)
    @test !has_active_cells_map(extended_grid)

    u = XFaceField(grid)
    set!(u, (x, y, z) -> z + 1e-5x)
    surface_u = view(u, :, :, Nz)

    extended_u = extended_external_variable(surface_u, extended_grid)
    Hx, Hy, _ = halo
    @test size(parent(extended_u)) == (Nx + 2Hx, Ny + 2Hy, 1)
    @test sizeof(parent(extended_u)) == (Nx + 2Hx) * (Ny + 2Hy) * sizeof(eltype(grid))
    @test interior(extended_u) == interior(surface_u)

    # Refreshing from the source keeps the copy in step, and its halo can be filled
    set!(u, (x, y, z) -> 2z - 1e-5x)
    refresh_and_fill_external_field!(extended_u, surface_u)
    @test interior(extended_u) == interior(surface_u)
end
