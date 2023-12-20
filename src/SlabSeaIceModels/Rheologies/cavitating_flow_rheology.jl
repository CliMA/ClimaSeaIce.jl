using Oceananigans.TimeSteppers: store_field_tendencies!

struct CavitatingFlowRheology{P, FT} <: AbstractExplicitRheology
    p  :: P
    P★ :: FT
    C  :: FT
    substeps :: Int
end

"""
    function CavitatingFlowRheology(grid::AbstractGrid; P★ = 1, C = 20, substeps = 100)

construct a `CavitatingFlowRheology` type, where ``∇ ⋅ σ = ∇ ⋅ p`` with 

- ``p = Pₚ`` if ``∇ ⋅ u < 0``
- ``p = 0`` otherwise

and ``Pₚ`` is the ice strength parameterized as ``P★ h exp( - C ⋅ ( 1 - ℵ ))``

Arguments
=========

- `grid`: the `SlabSeaIceModel` grid

Keyword Arguments
=================

- `P★`: ice strength multiplier (in Nm²)
- `C`: exponent coefficient
- `substeps`: number of substeps

"""
function CavitatingFlowRheology(grid::AbstractGrid; P★ = 1, C = 20, substeps = 100)
    p = CenterField(grid)
    FT = eltype(grid)
    return CavitatingFlowRheology(p, convert(FT, P★), convert(FT, C), substeps)
end

function compute_stresses!(model, rheology::CavitatingFlowRheology, Δτ) 

    grid = model.grid
    arch = architecture(grid)
    p = rheology.p
    C = rheology.C
    P★ = rheology.P★

    h = model.thickness
    ℵ = model.concentration

    u, v = model.velocities

    launch!(arch, grid, :xyz, _compute_cavitating_pressure!, p, grid, u, v, h, ℵ, P★, C)

    return nothing
end
   
@kernel function _compute_cavitating_pressure!(p, grid, u, v, h, ℵ, P★, C)
    i, j = @index(Global, NTuple)

    δ = div_xyᶜᶜᶜ(i, j, 1, grid, u, v)

    @inbounds p[i, j, 1] = ifelse(δ < 0, P★ * h[i, j, 1] * exp(- C * (1 - ℵ[i, j, 1])), 0)
end

@inline x_internal_stress_divergence(i, j, grid, ::CavitatingFlowRheology) = ∂xᶠᶜᶜ(i, j, 1, grid, rheology.p)
@inline y_internal_stress_divergence(i, j, grid, ::CavitatingFlowRheology) = ∂yᶜᶠᶜ(i, j, 1, grid, rheology.p)