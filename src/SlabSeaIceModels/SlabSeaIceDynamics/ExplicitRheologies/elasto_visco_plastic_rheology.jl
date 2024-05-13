using Oceananigans.TimeSteppers: store_field_tendencies!
using Oceananigans.Operators

import Oceananigans.BoundaryConditions: fill_halo_regions!

## The equations are solved in an iterative form following the EVP rheology of
## Kimmritz et al (2016) (https://www.sciencedirect.com/science/article/pii/S1463500317300690)
#
# Where:
# σᵢⱼ(u) = 2η ϵ̇ᵢⱼ + [(ζ - η) * (ϵ̇₁₁ + ϵ̇₂₂) - P / 2] δᵢⱼ
# uᵖ⁺¹ - uᵖ = β⁻¹ * (Δt / mᵢ * (∇ ⋅ σᵖ⁺¹ + fk̂ × uᵖ + τₐ + τₒ) + uⁿ - uᵖ)
#
struct ElastoViscoPlasticRheology{S1, S2, S3, U, V, P, FT, A} <: AbstractExplicitRheology
    σ₁₁   :: S1 # internal stress xx
    σ₂₂   :: S2 # internal stress yy
    σ₁₂   :: S3 # internal stress xy
    uⁿ    :: U  # ice u-velocity at time step n
    vⁿ    :: V  # ice v-velocity at time step n
    ice_strength :: P 
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    yield_curve_eccentricity :: FT # elliptic yield curve eccentricity
    ocean_ice_drag_coefficient :: FT
    Δ_min :: FT # minimum plastic parameter (transitions to viscous behaviour)
    substepping_coefficient :: A
    substeps :: Int
end

"""
    ElastoViscoPlasticRheology(grid::AbstractGrid; 
                               ice_compressive_strength = 27500, 
                               ice_compaction_hardening = 20, 
                               yield_curve_eccentricity = 2, 
                               ocean_ice_drag_coefficient = 5.5e-3,
                               Δ_min = 2e-9,
                               substepping_coefficient = ModifiedEVPSteppingCoefficients(grid),
                               substeps = 1000)

Constructs an `ElastoViscoPlasticRheology` object representing a "modified" elasto-visco-plastic
rheology for slab sea ice dynamics that follows the implementation of kimmritz et al (2016).
The `ice_strength` is parameterized as ``P★ h exp( - C ⋅ ( 1 - ℵ ))`` 

Arguments
=========
    
- `grid`: the `SlabSeaIceModel` grid

Keyword Arguments
=================
    
- `ice_compressive_strength`: parameter expressing compressive strength (in Nm²), default `27500`.
- `ice_compaction_hardening`: exponent coefficient for compaction hardening, default `20`.
- `yield_curve_eccentricity`: eccentricity of the elliptic yield curve, default `2`.
- `ocean_ice_drag_coefficient`: coefficient for the ocean - ice drag, default `5.5e-3`.
- `Δ_min`: Minimum value for the visco-plastic parameter. Limits the maximum viscosity of the ice, 
           transitioning the ice from a plastic to a viscous behaviour. Default value is `1e-10`.
- `substepping_coefficient`: Coefficient for substepping momentum (β) and internal stresses (α) (depends on the particular EVP formulation).
              Default value is `ModifiedEVPSteppingCoefficients(grid)`.
- `substeps`: Number of substeps for the visco-plastic calculation. Default value is `100`.
              Note that we here assume that β (the modified EVP parameter that substeps velocity)
              is equal to the number of substeps.
"""
function ElastoViscoPlasticRheology(grid::AbstractGrid; 
                                    ice_compressive_strength = 27500, 
                                    ice_compaction_hardening = 20, 
                                    yield_curve_eccentricity = 2, 
                                    ocean_ice_drag_coefficient = 5.5e-3,
                                    Δ_min = 2e-9,
                                    substepping_coefficient = ModifiedEVPSteppingCoefficients(grid),
                                    substeps = 1000)
    σ₁₁ = CenterField(grid)
    σ₂₂ = CenterField(grid)
    σ₁₂ = Field{Face, Face, Center}(grid)
    uⁿ = XFaceField(grid) 
    vⁿ = YFaceField(grid) 

    P  = CenterField(grid)
    FT = eltype(grid)
    return ElastoViscoPlasticRheology(σ₁₁, σ₂₂, σ₁₂, uⁿ, vⁿ, P, 
                                      convert(FT, ice_compressive_strength), 
                                      convert(FT, ice_compaction_hardening), 
                                      convert(FT, yield_curve_eccentricity),
                                      convert(FT, ocean_ice_drag_coefficient),
                                      convert(FT, Δ_min),
                                      substepping_coefficient,
                                      substeps)
end

# Extend the `adapt_structure` function for the ElastoViscoPlasticRheology
Adapt.adapt_structure(to, r::ElastoViscoPlasticRheology) = 
    ElastoViscoPlasticRheology(Adapt.adapt(to, r.σ₁₁),
                               Adapt.adapt(to, r.σ₂₂),
                               Adapt.adapt(to, r.σ₁₂),
                               Adapt.adapt(to, r.uⁿ),
                               Adapt.adapt(to, r.vⁿ),
                               Adapt.adapt(to, r.ice_strength),
                               Adapt.adapt(to, r.ice_compressive_strength),
                               Adapt.adapt(to, r.ice_compaction_hardening),
                               Adapt.adapt(to, r.yield_curve_eccentricity),
                               Adapt.adapt(to, r.ocean_ice_drag_coefficient),
                               Adapt.adapt(to, r.Δ_min),
                               Adapt.adapt(to, r.substepping_coefficient),
                               r.substeps)

# Extend `fill_halo_regions!` for the ElastoViscoPlasticRheology
function fill_halo_regions!(rheology::ElastoViscoPlasticRheology)
    fill_halo_regions!(rheology.σ₁₁)
    fill_halo_regions!(rheology.σ₂₂)
    fill_halo_regions!(rheology.σ₁₂)

    return nothing
end

"""
    initialize_substepping!(model, rheology::ElastoViscoPlasticRheology)

Initialize the substepping for the elasto-visco-plastic rheology.
In this step we save down the velocities at the previous time step and
calculate the ice strength given the ice mass (thickness and concentration).
"""
function initialize_substepping!(model, rheology::ElastoViscoPlasticRheology)
    uⁿ = model.velocities.u
    vⁿ = model.velocities.v
    h  = model.ice_thickness
    ℵ  = model.concentration

    launch!(architecture(model.grid), model.grid, :xy, _store_initial_velocities!, rheology, uⁿ, vⁿ)
    launch!(architecture(model.grid), model.grid, :xy, _compute_ice_strength!,     rheology, h, ℵ)

    fill_halo_regions!(rheology.uⁿ)
    fill_halo_regions!(rheology.vⁿ)
    fill_halo_regions!(rheology.ice_strength)

    return nothing
end

# We need initial velocities for the momentum update step
@kernel function _store_initial_velocities!(rheology, uⁿ, vⁿ)
    i, j = @index(Global, NTuple)

    u, v = rheology.uⁿ, rheology.vⁿ
    @inbounds u[i, j, 1] = uⁿ[i, j, 1]
    @inbounds v[i, j, 1] = vⁿ[i, j, 1]
end

# Specific compute stresses for the EVP rheology
function compute_stresses!(model, rheology::ElastoViscoPlasticRheology, Δt) 

    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ℵ  = model.concentration
    ρᵢ = model.ice_density
    
    u, v = model.velocities
    launch!(arch, grid, :xyz, _compute_evp_stresses!, rheology, grid, u, v, h, ℵ, ρᵢ, Δt)

    return nothing
end

# Compute the elasto-visco-plastic stresses for a slab sea ice model.
# The function updates the internal stress variables `σ₁₁`, `σ₂₂`, and `σ₁₂` in the `rheology` object
# following the mEVP formulation of Kimmritz et al (2016).
# This is the `meat` of the formulation.
@kernel function _compute_evp_stresses!(rheology, grid, u, v, h, ℵ, ρᵢ, Δt)
    i, j = @index(Global, NTuple)

    P   = rheology.ice_strength
    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.Δ_min

    # Extract internal stresses
    σ₁₁ = rheology.σ₁₁
    σ₂₂ = rheology.σ₂₂
    σ₁₂ = rheology.σ₁₂

    # Strain rates
    ϵ̇₁₁ =  ∂xᶜᶜᶜ(i, j, 1, grid, u)
    ϵ̇₁₂ = (∂xᶠᶠᶜ(i, j, 1, grid, v) + ∂yᶠᶠᶜ(i, j, 1, grid, u)) / 2
    ϵ̇₂₂ =  ∂yᶜᶜᶜ(i, j, 1, grid, v)

    # Center - Center variables:
    ϵ̇₁₂ᶜᶜᶜ = (ℑxyᶜᶜᶜ(i, j, 1, grid, ∂xᶠᶠᶜ, v) + ℑxyᶜᶜᶜ(i, j, 1, grid, ∂yᶠᶠᶜ, u)) / 2

    # Ice divergence 
    δ = ϵ̇₁₁ + ϵ̇₂₂

    # Ice shear (at Centers)
    s = sqrt((ϵ̇₁₁ - ϵ̇₂₂)^2 + 4ϵ̇₁₂ᶜᶜᶜ^2)

    # Visco - Plastic parameter 
    # if Δ is very small we assume a linear viscous response
    # adding a minimum Δ_min (at Centers)
    Δᶜᶜᶜ = sqrt(δ^2 + s^2 * e⁻²) + Δm

    # Face - Face variables
    ϵ̇₁₁ᶠᶠᶜ = ℑxyᶠᶠᶜ(i, j, 1, grid, ∂xᶜᶜᶜ, u)
    ϵ̇₂₂ᶠᶠᶜ = ℑxyᶠᶠᶜ(i, j, 1, grid, ∂yᶜᶜᶜ, v)

    # Ice divergence
    δᶠᶠᶜ = ϵ̇₁₁ᶠᶠᶜ + ϵ̇₂₂ᶠᶠᶜ

    # Ice shear
    sᶠᶠᶜ = sqrt((ϵ̇₁₁ᶠᶠᶜ - ϵ̇₂₂ᶠᶠᶜ)^2 + 4ϵ̇₁₂^2)

    # Visco - Plastic parameter 
    # if Δ is very small we assume a linear viscous response
    # adding a minimum Δ_min (at Faces)
    Δᶠᶠᶜ = sqrt(δᶠᶠᶜ^2 + sᶠᶠᶜ^2 * e⁻²) + Δm

    # ice strength calculation 
    # Note: can we interpolate P on faces or do we need to compute it on faces?
    Pᶜᶜᶜ = @inbounds P[i, j, 1]
    Pᶠᶠᶜ = ℑxyᶠᶠᶜ(i, j, 1, grid, P)

    # ζ: Bulk viscosity (viscosity which responds to compression) 
    # η: Shear viscosity (viscosity which responds to shear)
    ζᶜᶜᶜ = Pᶜᶜᶜ / (2Δᶜᶜᶜ)
    ηᶜᶜᶜ = ζᶜᶜᶜ * e⁻²

    ζᶠᶠᶜ = Pᶠᶠᶜ / (2Δᶠᶠᶜ)
    ηᶠᶠᶜ = ζᶠᶠᶜ * e⁻²

    # Replacement pressure
    Pᵣ = Pᶜᶜᶜ * Δᶜᶜᶜ / (Δᶜᶜᶜ + Δm)

    # σ(uᵖ): the tangential stress depends only shear viscosity 
    # while the compressive stresses depend on the bulk viscosity and the ice strength
    σ₁₁ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₁₁ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2) 
    σ₂₂ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₂₂ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2)
    σ₁₂ᵖ⁺¹ = 2 * ηᶠᶠᶜ * ϵ̇₁₂

    mᵢ = @inbounds h[i, j, 1] * ℵ[i, j, 1] * ρᵢ

    c = rheology.substepping_coefficient

    # Update coefficients for substepping if we are using dynamic substepping
    # with spatially varying coefficients such as in Kimmritz et al (2016)
    update_stepping_coefficients!(i, j, 1, grid, c, ζᶜᶜᶜ, mᵢ, Δt)

    # Coefficient for substepping internal stress
    α = get_stepping_coefficients(i, j, 1, grid, rheology, c)

    @inbounds σ₁₁[i, j, 1] += ifelse(mᵢ > 0, (σ₁₁ᵖ⁺¹ - σ₁₁[i, j, 1]) / α, 0)
    @inbounds σ₂₂[i, j, 1] += ifelse(mᵢ > 0, (σ₂₂ᵖ⁺¹ - σ₂₂[i, j, 1]) / α, 0)
    @inbounds σ₁₂[i, j, 1] += ifelse(mᵢ > 0, (σ₁₂ᵖ⁺¹ - σ₁₂[i, j, 1]) / α, 0)
end

#####
##### Internal stress divergence for the EVP model
#####

@inline function x_internal_stress_divergence(i, j, grid, r::ElastoViscoPlasticRheology) 
    ∂xσ₁₁ = δxᶠᶜᶜ(i, j, 1, grid, Ax_qᶜᶜᶜ, r.σ₁₁)
    ∂yσ₁₂ = δyᶠᶜᶜ(i, j, 1, grid, Ay_qᶠᶠᶜ, r.σ₁₂)

    return (∂xσ₁₁ + ∂yσ₁₂) / Vᶠᶜᶜ(i, j, 1, grid)
end

@inline function y_internal_stress_divergence(i, j, grid, r::ElastoViscoPlasticRheology) 
    ∂xσ₁₂ = δxᶜᶠᶜ(i, j, 1, grid, Ax_qᶠᶠᶜ, r.σ₁₂)
    ∂yσ₂₂ = δyᶜᶠᶜ(i, j, 1, grid, Ay_qᶜᶜᶜ, r.σ₂₂)

    return (∂xσ₁₂ + ∂yσ₂₂) / Vᶜᶠᶜ(i, j, 1, grid)
end
