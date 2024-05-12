using Oceananigans.TimeSteppers: store_field_tendencies!

struct CavitatingFlowRheology{P, I, FT} <: AbstractExplicitRheology
    pressure :: P
    ice_strength :: I
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT
    substeps :: Int
end

"""
    function CavitatingFlowRheology(grid::AbstractGrid; P★ = 1, C = 20, substeps = 100)

construct a `CavitatingFlowRheology` type, where ``∇ ⋅ σ = ∇p`` with 

- ``p = Pₚ`` if ``∇ ⋅ u < 0``
- ``p = 0`` otherwise

and ``Pₚ`` is the ice strength parameterized as ``P★ h exp( - C ⋅ ( 1 - ℵ ))``

Arguments
=========

- `grid`: the `SlabSeaIceModel` grid

Keyword Arguments
=================

- `P★`: parameter expressing compressive strength (in Nm²)
- `C`: exponent coefficient for compaction hardening
- `substeps`: number of substeps

"""
function CavitatingFlowRheology(grid::AbstractGrid; 
                                ice_compressive_strength = 27500, 
                                ice_compaction_hardening = 20,
                                substeps = 100)
    p = CenterField(grid)
    P = CenterField(grid)
    FT = eltype(grid)
    return CavitatingFlowRheology(p, P, 
                                  convert(FT, ice_compressive_strength), 
                                  convert(FT, ice_compaction_hardening), 
                                  substeps)
end

function initialize_substepping!(model, rheology::CavitatingFlowRheology)
    h = model.ice_thickness
    ℵ = model.concentration

    launch!(architecture(model.grid), model.grid, :xy, _compute_ice_strength!, rheology, h, ℵ)
    
    fill_halo_regions!(rheology.ice_strength)

    return nothing
end

@kernel function _compute_ice_strength!(rheology, h, ℵ)
    i, j = @index(Global, NTuple)
    
    P    = rheology.ice_strength
    P★   = rheology.ice_compressive_strength
    C    = rheology.ice_compaction_hardening

    @inbounds P[i, j, 1] = @inbounds P★ * h[i, j, 1] * exp(- C * (1 - ℵ[i, j, 1])) 
end

function compute_stresses!(model, rheology::CavitatingFlowRheology, Δt) 

    grid = model.grid
    arch = architecture(grid)
    p = rheology.pressure
    P = rheology.ice_strength

    h = model.thickness
    ℵ = model.concentration

    u, v = model.velocities

    launch!(arch, grid, :xyz, _compute_cavitating_pressure!, p, grid, u, v, h, ℵ, P)

    return nothing
end
   
@kernel function _compute_cavitating_pressure!(p, grid, u, v, h, ℵ, P)
    i, j = @index(Global, NTuple)

    δ = div_xyᶜᶜᶜ(i, j, 1, grid, u, v)

    @inbounds p[i, j, 1] = ifelse(δ <= 0, P[i, j, k], 0)
end

@inline x_internal_stress_divergence(i, j, grid, r::CavitatingFlowRheology) = - ∂xᶠᶜᶜ(i, j, 1, grid, r.p)
@inline y_internal_stress_divergence(i, j, grid, r::CavitatingFlowRheology) = - ∂yᶜᶠᶜ(i, j, 1, grid, r.p)