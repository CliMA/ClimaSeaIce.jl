using Oceananigans.Architectures: architecture, on_architecture
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

SeaIceColumnDiscretization(reference_faces) =
    SeaIceColumnDiscretization(reference_faces, reference_faces, nothing, nothing, nothing, nothing, nothing, nothing)

function Grids.validate_dimension_specification(T, ξ::SeaIceColumnDiscretization, dir, N, FT)
    cᵃᵃᶠ = Grids.validate_dimension_specification(T, ξ.cᵃᵃᶠ, dir, N, FT)
    cᵃᵃᶜ = Grids.validate_dimension_specification(T, ξ.cᵃᵃᶜ, dir, N, FT)
    return SeaIceColumnDiscretization(cᵃᵃᶠ, cᵃᵃᶜ, ξ.Δᵃᵃᶠ, ξ.Δᵃᵃᶜ, ξ.hsⁿ, ξ.hbⁿ, ξ.hs⁻, ξ.hb⁻)
end

function Adapt.adapt_structure(to, coordinate::SeaIceColumnDiscretization)
    return SeaIceColumnDiscretization(Adapt.adapt(to, coordinate.cᵃᵃᶠ),
                                      Adapt.adapt(to, coordinate.cᵃᵃᶜ),
                                      Adapt.adapt(to, coordinate.Δᵃᵃᶠ),
                                      Adapt.adapt(to, coordinate.Δᵃᵃᶜ),
                                      Adapt.adapt(to, coordinate.hsⁿ),
                                      Adapt.adapt(to, coordinate.hbⁿ),
                                      Adapt.adapt(to, coordinate.hs⁻),
                                      Adapt.adapt(to, coordinate.hb⁻))
end

# Materialize the reference layers and allocate the two-interface height fields during grid construction.
function Grids.generate_coordinate(FT, topology, size, halo, coordinate::SeaIceColumnDiscretization, coordinate_name, dim::Int, arch)
    dim == 3 || throw(ArgumentError("SeaIceColumnDiscretization is supported only in the third dimension (z)"))
    coordinate_name == :z ||
        throw(ArgumentError("SeaIceColumnDiscretization is supported only for the z-coordinate"))

    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
    reference_faces = coordinate.cᵃᵃᶠ

    Lr, rᵃᵃᶠ, rᵃᵃᶜ, Δrᵃᵃᶠ, Δrᵃᵃᶜ = Grids.generate_coordinate(FT, topology[3](), Nz, Hz, reference_faces, :r, arch)

    args = (FT, arch, (Center, Center, Nothing), topology, (Nx, Ny, Nz), (Hx, Hy, Hz))
    hsⁿ = new_data(args...)
    hbⁿ = new_data(args...)
    hs⁻ = new_data(args...)
    hb⁻ = new_data(args...)

    # Default the interfaces to a resting column of the reference height (base at 0, surface at Lr), so a fresh
    # grid behaves like a static one until `initialize_column_interfaces!`/the coupled step move the interfaces.
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

@inline function column_height(grid, i, j)
    z = sea_ice_discretization(grid)
    return @inbounds z.hsⁿ[i, j, 1] - z.hbⁿ[i, j, 1]
end

@inline function previous_column_height(grid, i, j)
    z = sea_ice_discretization(grid)
    return @inbounds z.hs⁻[i, j, 1] - z.hb⁻[i, j, 1]
end

# The vertical metric is the column height divided by the reference height (the z-star convention). Defining σⁿ/σ⁻
# lets the standard Oceananigans operators — `Δz = Δr σ`, `znode = r σ + η`, the metric ratio `σ⁻/σⁿ` — reuse it,
# and Oceananigans forwards σⁿ/σ⁻ from an immersed grid to its underlying grid automatically.
@inline Operators.σⁿ(i, j, k, grid::SeaIceColumnUnderlyingGrid, ℓx, ℓy, ℓz) = column_height(grid, i, j) / grid.Lz
@inline Operators.σ⁻(i, j, k, grid::SeaIceColumnUnderlyingGrid, ℓx, ℓy, ℓz) = previous_column_height(grid, i, j) / grid.Lz

@inline Operators.Δzᶜᶜᶜ(i, j, k, grid::SeaIceColumnGrid) = Operators.Δrᶜᶜᶜ(i, j, k, grid) * Operators.σⁿ(i, j, k, grid, Center(), Center(), Center())
@inline Operators.Δzᶜᶜᶠ(i, j, k, grid::SeaIceColumnGrid) = Operators.Δrᶜᶜᶠ(i, j, k, grid) * Operators.σⁿ(i, j, k, grid, Center(), Center(), Face())

@inline Grids.znode(i, j, k, grid::SeaIceColumnUnderlyingGrid, ::Center, ::Center, ℓz) = rnode(i, j, k, grid, Center(), Center(), ℓz) * Operators.σⁿ(i, j, k, grid, Center(), Center(), ℓz) + @inbounds grid.z.hbⁿ[i, j, 1]

# Place every column at rest: surface at z = 0, base at z = -ice_thickness, previous heights matching.
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
    launch!(architecture(grid), grid, :xy, _initialize_column_interfaces!,
            sea_ice_discretization(grid), ice_thickness)
    return nothing
end

# Roll current interface heights into the previous slot, then move the base by `basal_growth` (positive grows the
# ice downward, negative melts it upward) and the surface by `surface_growth` (positive for snow-ice, negative for
# surface melt). The retained previous heights make this step's motion visible to the swept-face displacement.
@kernel function _advance_column_interfaces!(z, basal_growth, surface_growth)
    i, j = @index(Global, NTuple)
    @inbounds begin
        z.hs⁻[i, j, 1] = z.hsⁿ[i, j, 1]
        z.hb⁻[i, j, 1] = z.hbⁿ[i, j, 1]
        z.hsⁿ[i, j, 1] += surface_growth[i, j, 1]
        z.hbⁿ[i, j, 1] -= basal_growth[i, j, 1]
    end
end

function advance_column_interfaces!(grid, basal_growth, surface_growth)
    launch!(architecture(grid), grid, :xy, _advance_column_interfaces!, sea_ice_discretization(grid), basal_growth, surface_growth)
    return nothing
end
