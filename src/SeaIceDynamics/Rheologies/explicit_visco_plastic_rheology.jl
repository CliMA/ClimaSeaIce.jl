using Oceananigans.TimeSteppers: store_field_tendencies!
using Oceananigans.Operators
using Oceananigans.Grids: AbstractGrid

## The equations are solved in an iterative form following the EVP rheology of
## Kimmritz et al (2016) (https://www.sciencedirect.com/science/article/pii/S1463500317300690)
#
# Where:
# σᵢⱼ(u) = 2η ϵ̇ᵢⱼ + [(ζ - η) * (ϵ̇₁₁ + ϵ̇₂₂) - P / 2] δᵢⱼ
#
struct ExplicitViscoPlasticRheology{S1, S2, S3, U, P, FT} <: AbstractRheology
    σ₁₁   :: S1 # internal stress xx
    σ₂₂   :: S2 # internal stress yy
    σ₁₂   :: S3 # internal stress xy
    previous_velocities :: U # ice velocities at time step n
    ice_strength :: P # field containing the precomputed ice strength
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    yield_curve_eccentricity :: FT # elliptic yield curve eccentricity
    Δ_min :: FT # minimum plastic parameter (transitions to viscous behaviour)
end

"""
    ExplicitViscoPlasticRheology(grid::AbstractGrid; 
                               ice_compressive_strength = 27500, 
                               ice_compaction_hardening = 20, 
                               yield_curve_eccentricity = 2, 
                               ocean_ice_drag_coefficient = 5.5e-3,
                               Δ_min = 2e-9,
                               substepping_coefficient = DynamicSteppingCoefficient(grid),
                               substeps = 1000)

Constructs an `ExplicitViscoPlasticRheology` object representing a "modified" elasto-visco-plastic
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
- `Δ_min`: Minimum value for the visco-plastic parameter. Limits the maximum viscosity of the ice, 
           transitioning the ice from a plastic to a viscous behaviour. Default value is `1e-10`.
"""
function ExplicitViscoPlasticRheology(grid::AbstractGrid; 
                                      ice_compressive_strength = 27500, 
                                      ice_compaction_hardening = 20, 
                                      yield_curve_eccentricity = 2, 
                                      Δ_min = 2e-9)
    σ₁₁ = CenterField(grid)
    σ₂₂ = CenterField(grid)
    σ₁₂ = Field{Face, Face, Center}(grid)

    uⁿ = XFaceField(grid) 
    vⁿ = YFaceField(grid) 

    previous_velocities = (u = uⁿ, v = vⁿ)

    P  = CenterField(grid)
    FT = eltype(grid)
    return ExplicitViscoPlasticRheology(σ₁₁, σ₂₂, σ₁₂, P, previous_velocities,
                                      convert(FT, ice_compressive_strength), 
                                      convert(FT, ice_compaction_hardening), 
                                      convert(FT, yield_curve_eccentricity),
                                      convert(FT, Δ_min))
end

# Extend the `adapt_structure` function for the ExplicitViscoPlasticRheology
Adapt.adapt_structure(to, r::ExplicitViscoPlasticRheology) = 
    ExplicitViscoPlasticRheology(Adapt.adapt(to, r.σ₁₁),
                                 Adapt.adapt(to, r.σ₂₂),
                                 Adapt.adapt(to, r.σ₁₂),
                                 Adapt.adapt(to, r.previous_velocities),
                                 Adapt.adapt(to, r.ice_strength),
                                 Adapt.adapt(to, r.ice_compressive_strength),
                                 Adapt.adapt(to, r.ice_compaction_hardening),
                                 Adapt.adapt(to, r.yield_curve_eccentricity),
                                 Adapt.adapt(to, r.Δ_min))

"""
    initialize_rheology!(model, rheology::ExplicitViscoPlasticRheology)

Initialize the elasto-visco-plastic rheology.
In this step we calculate the ice strength given the ice mass (thickness and concentration).
"""
function initialize_rheology!(model, rheology::ExplicitViscoPlasticRheology)
    h = model.ice_thickness
    ℵ = model.ice_concentration

    P  = rheology.ice_strength
    P★ = rheology.ice_compressive_strength
    C  = rheology.ice_compaction_hardening
    
    u, v   = model.velocities
    uⁿ, vⁿ = rheology.previous_velocities

    # compute on the whole grid including halos
    parameters = KernelParameters(size(P.data), P.data.offsets)
    launch!(architecture(model.grid), model.grid, parameters, _initialize_evp_rhology!, P, P★, C, h, ℵ, uⁿ, vⁿ, u, v)
    
    return nothing
end

# Extend `fill_halo_regions!` for the AbstractRheology
function fill_halo_regions!(rheology::ExplicitViscoPlasticRheology, args...)
    fill_halo_regions!(rheology.σ₁₁, args...)
    fill_halo_regions!(rheology.σ₂₂, args...)
    fill_halo_regions!(rheology.σ₁₂, args...)

    return nothing
end

# Extend `mask_immersed_field!` for the AbstractRheology
function mask_immersed_field!(rheology::ExplicitViscoPlasticRheology)
    mask_immersed_field!(rheology.σ₁₁)
    mask_immersed_field!(rheology.σ₂₂)
    mask_immersed_field!(rheology.σ₁₂)

    return nothing
end

@kernel function _initialize_evp_rhology!(P, P★, C, h, ℵ, uⁿ, vⁿ, u, v)
    i, j = @index(Global, NTuple)    
    @inbounds  P[i, j, 1] = ice_strength(P★, C, h, ℵ)
    @inbounds uⁿ[i, j, 1] = u[i, j, 1]
    @inbounds vⁿ[i, j, 1] = v[i, j, 1]
end

# The parameterization for an `ExplicitViscoPlasticRheology`
@inline ice_strength(P★, C, h, ℵ) = P★ * h[i, j, 1] * exp(- C * (1 - ℵ[i, j, 1])) 

# Specific compute stresses for the EVP rheology
function compute_stresses!(model, solver, rheology::ExplicitViscoPlasticRheology, Δt) 

    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ℵ  = model.ice_concentration
    ρᵢ = model.ice_density

    substeps = solver.substeps
    stepping_coefficient = solver.substepping_coefficient
    
    u, v = model.velocities
    launch!(arch, grid, :xyz, _compute_evp_stresses!, rheology, grid, u, v, h, ℵ, ρᵢ, Δt, substeps, stepping_coefficient)

    return nothing
end

# Compute the visco-plastic stresses for a slab sea ice model.
# The function updates the internal stress variables `σ₁₁`, `σ₂₂`, and `σ₁₂` in the `rheology` object
# following the mEVP formulation of Kimmritz et al (2016).
# This is the `meat` of the formulation.
@kernel function _compute_evp_stresses!(rheology, grid, u, v, h, ℵ, ρᵢ, Δt, substeps, stepping_coefficient)
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
    ϵ̇₁₂ᶜᶜᶜ = (ℑxyᴮᶜᶜᶜ(i, j, 1, grid, ∂xᶠᶠᶜ, v) + ℑxyᴮᶜᶜᶜ(i, j, 1, grid, ∂yᶠᶠᶜ, u)) / 2

    # Ice divergence 
    δ = ϵ̇₁₁ + ϵ̇₂₂

    # Ice shear (at Centers)
    s = sqrt((ϵ̇₁₁ - ϵ̇₂₂)^2 + 4ϵ̇₁₂ᶜᶜᶜ^2)

    # Visco - Plastic parameter 
    # if Δ is very small we assume a linear viscous response
    # adding a minimum Δ_min (at Centers)
    Δᶜᶜᶜ = sqrt(δ^2 + s^2 * e⁻²) + Δm

    # Face - Face variables
    ϵ̇₁₁ᶠᶠᶜ = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, ∂xᶜᶜᶜ, u)
    ϵ̇₂₂ᶠᶠᶜ = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, ∂yᶜᶜᶜ, v)

    # Ice divergence
    δᶠᶠᶜ = ϵ̇₁₁ᶠᶠᶜ + ϵ̇₂₂ᶠᶠᶜ

    # Ice shear
    sᶠᶠᶜ = sqrt((ϵ̇₁₁ᶠᶠᶜ - ϵ̇₂₂ᶠᶠᶜ)^2 + 4ϵ̇₁₂^2)

    # Visco - Plastic parameter 
    # if Δ is very small we assume a linear viscous response
    # adding a minimum Δ_min (at Faces)
    Δᶠᶠᶜ = sqrt(δᶠᶠᶜ^2 + sᶠᶠᶜ^2 * e⁻²) + Δm

    # Ice strength calculation 
    # Note: can we interpolate P on faces or do we need to compute it on faces?
    Pᶜᶜᶜ = @inbounds P[i, j, 1]
    Pᶠᶠᶜ = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, P)

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

    c = stepping_coefficient

    # Update coefficients for substepping if we are using dynamic substepping
    # with spatially varying coefficients such as in Kimmritz et al (2016)
    update_stepping_coefficients!(i, j, 1, grid, c, ζᶜᶜᶜ, mᵢ, Δt)

    # Coefficient for substepping internal stress
    α = get_stepping_coefficients(i, j, 1, grid, substeps, c)

    @inbounds σ₁₁[i, j, 1] += ifelse(mᵢ > 0, (σ₁₁ᵖ⁺¹ - σ₁₁[i, j, 1]) / α, zero(grid))
    @inbounds σ₂₂[i, j, 1] += ifelse(mᵢ > 0, (σ₂₂ᵖ⁺¹ - σ₂₂[i, j, 1]) / α, zero(grid))
    @inbounds σ₁₂[i, j, 1] += ifelse(mᵢ > 0, (σ₁₂ᵖ⁺¹ - σ₁₂[i, j, 1]) / α, zero(grid))
end

#####
##### Internal stress divergence for the EVP model
#####

@inline function x_internal_stress_divergence(i, j, k, grid, r::ExplicitViscoPlasticRheology) 
    ∂xσ₁₁ = δxᶠᶜᶜ(i, j, k, grid, Ax_qᶜᶜᶜ, r.σ₁₁)
    ∂yσ₁₂ = δyᶠᶜᶜ(i, j, k, grid, Ay_qᶠᶠᶜ, r.σ₁₂)

    return (∂xσ₁₁ + ∂yσ₁₂) / Vᶠᶜᶜ(i, j, k, grid)
end

@inline function y_internal_stress_divergence(i, j, k, grid, r::ExplicitViscoPlasticRheology) 
    ∂xσ₁₂ = δxᶜᶠᶜ(i, j, k, grid, Ax_qᶠᶠᶜ, r.σ₁₂)
    ∂yσ₂₂ = δyᶜᶠᶜ(i, j, k, grid, Ay_qᶜᶜᶜ, r.σ₂₂)

    return (∂xσ₁₂ + ∂yσ₂₂) / Vᶜᶠᶜ(i, j, k, grid)
end

# To help convergence to the right velocities
@inline rheology_specific_numerical_terms_x(i, j, k, grid, r::ExplicitViscoPlasticRheology, uᵢ) = r.uⁿ[i, j, k] - uᵢ[i, j, k]
@inline rheology_specific_numerical_terms_y(i, j, k, grid, r::ExplicitViscoPlasticRheology, vᵢ) = r.vⁿ[i, j, k] - vᵢ[i, j, k]