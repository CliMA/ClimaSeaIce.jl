using Oceananigans.Architectures: architecture
using Oceananigans.Fields: Field, Center, Face
using Oceananigans.Grids: Grids, AbstractVerticalCoordinate, AbstractUnderlyingGrid, Bounded, rnode, new_data
using Oceananigans.Operators: Operators
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index

"""
    SeaIceColumnDiscretization

Vertical coordinate for a sea-ice column with two independently moving interfaces. The reference layers partition
`[r_bottom, r_top]` and the physical height of reference node `r` in column `(i, j)` is

    z(r) = hbⁿ + (r - r_bottom) / (r_top - r_bottom) * (hsⁿ - hbⁿ),

so the column thickness is `hsⁿ - hbⁿ`. The base interface moves with basal congelation and melt;
the surface interface moves with surface melt and snow-ice formation. The previous-step interface heights are retained
so the per-face vertical displacement — and hence the conservative swept-face enthalpy flux — can distinguish growth or
melt at the base from that at the surface.

`cᵃᵃᶠ`/`cᵃᵃᶜ`/`Δᵃᵃᶠ`/`Δᵃᵃᶜ` are the reference (`r`) node/spacing arrays, named as Oceananigans expects so `Δrᶜᶜᶜ`
and `rnode` apply directly; the height fields are `nothing` until `generate_coordinate` allocates them during grid
construction.
"""
struct SeaIceColumnDiscretization{CF, CC, DF, DC, H} <: AbstractVerticalCoordinate
    cᵃᵃᶠ :: CF
    cᵃᵃᶜ :: CC
    Δᵃᵃᶠ :: DF
    Δᵃᵃᶜ :: DC
     hsⁿ :: H
     hbⁿ :: H
     hs⁻ :: H
     hb⁻ :: H
end

function SeaIceColumnDiscretization(reference_faces)
    return SeaIceColumnDiscretization(reference_faces, reference_faces, nothing, nothing, nothing, nothing, nothing, nothing)
end

function Grids.validate_dimension_specification(T, ξ::SeaIceColumnDiscretization, dir, N, FT)
    cᵃᵃᶠ = Grids.validate_dimension_specification(T, ξ.cᵃᵃᶠ, dir, N, FT)
    cᵃᵃᶜ = Grids.validate_dimension_specification(T, ξ.cᵃᵃᶜ, dir, N, FT)
    return SeaIceColumnDiscretization(cᵃᵃᶠ, cᵃᵃᶜ, ξ.Δᵃᵃᶠ, ξ.Δᵃᵃᶜ, ξ.hsⁿ, ξ.hbⁿ, ξ.hs⁻, ξ.hb⁻)
end

Adapt.@adapt_structure SeaIceColumnDiscretization

# Materialize the reference layers and allocate the two-interface height fields during grid construction.
function Grids.generate_coordinate(FT, topology, size, halo, coordinate::SeaIceColumnDiscretization, coordinate_name, dim::Int, arch)
    dim == 3 && coordinate_name == :z || throw(ArgumentError("SeaIceColumnDiscretization is supported only for the z-coordinate"))

    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
    reference_faces = coordinate.cᵃᵃᶠ

    Lr, rᵃᵃᶠ, rᵃᵃᶜ, Δrᵃᵃᶠ, Δrᵃᵃᶜ = Grids.generate_coordinate(FT, topology[3](), Nz, Hz, reference_faces, :r, arch)

    args = (FT, arch, (Center, Center, Nothing), topology, (Nx, Ny, Nz), (Hx, Hy, Hz))
    hsⁿ = new_data(args...)
    hbⁿ = new_data(args...)
    hs⁻ = new_data(args...)
    hb⁻ = new_data(args...)

    # A fresh grid is a resting column of the reference height: base at 0, surface at Lr.
    fill!(hbⁿ, zero(FT)); fill!(hb⁻, zero(FT))
    fill!(hsⁿ, convert(FT, Lr)); fill!(hs⁻, convert(FT, Lr))

    coordinate = SeaIceColumnDiscretization(rᵃᵃᶠ, rᵃᵃᶜ, Δrᵃᵃᶠ, Δrᵃᵃᶜ, hsⁿ, hbⁿ, hs⁻, hb⁻)

    return Lr, coordinate
end

const SeaIceColumnUnderlyingGrid = AbstractUnderlyingGrid{<:Any, <:Any, <:Any, <:Bounded, <:SeaIceColumnDiscretization}
const SeaIceColumnGrid = Union{SeaIceColumnUnderlyingGrid,
                               ImmersedBoundaryGrid{<:Any, <:Any, <:Any, <:Any, <:SeaIceColumnUnderlyingGrid}}

@inline sea_ice_discretization(grid::SeaIceColumnUnderlyingGrid) = grid.z
@inline sea_ice_discretization(grid::ImmersedBoundaryGrid) = grid.underlying_grid.z

@inline column_height(grid, i, j) = @inbounds sea_ice_discretization(grid).hsⁿ[i, j, 1] - sea_ice_discretization(grid).hbⁿ[i, j, 1]
@inline previous_column_height(grid, i, j) = @inbounds sea_ice_discretization(grid).hs⁻[i, j, 1] - sea_ice_discretization(grid).hb⁻[i, j, 1]

# z-star metric: σ = column height / reference height.
@inline Operators.σⁿ(i, j, k, grid::SeaIceColumnUnderlyingGrid, ℓx, ℓy, ℓz) = column_height(grid, i, j) / grid.Lz
@inline Operators.σ⁻(i, j, k, grid::SeaIceColumnUnderlyingGrid, ℓx, ℓy, ℓz) = previous_column_height(grid, i, j) / grid.Lz

@inline Operators.Δzᶜᶜᶜ(i, j, k, grid::SeaIceColumnGrid) = Operators.Δrᶜᶜᶜ(i, j, k, grid) * Operators.σⁿ(i, j, k, grid, Center(), Center(), Center())
@inline Operators.Δzᶜᶜᶠ(i, j, k, grid::SeaIceColumnGrid) = Operators.Δrᶜᶜᶠ(i, j, k, grid) * Operators.σⁿ(i, j, k, grid, Center(), Center(), Face())

@inline Grids.znode(i, j, k, grid::SeaIceColumnUnderlyingGrid, ::Center, ::Center, ℓz) = rnode(i, j, k, grid, Center(), Center(), ℓz) * Operators.σⁿ(i, j, k, grid, Center(), Center(), ℓz) + @inbounds grid.z.hbⁿ[i, j, 1]

# Columns at rest: surface at z = 0, base at z = -ice_thickness.
@kernel function _initialize_column_interfaces!(z, ice_thickness)
    i, j = @index(Global, NTuple)
    @inbounds begin
        h = ice_thickness[i, j, 1]
        z.hsⁿ[i, j, 1] = zero(h)
        z.hbⁿ[i, j, 1] = -h
        z.hs⁻[i, j, 1] = zero(h)
        z.hb⁻[i, j, 1] = -h
    end
end

function initialize_column_interfaces!(grid, ice_thickness)
    launch!(architecture(grid), grid, :xy, _initialize_column_interfaces!, sea_ice_discretization(grid), ice_thickness)
    return nothing
end
