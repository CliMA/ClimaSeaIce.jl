using Oceananigans.Operators
using Oceananigans.Grids: AbstractGrid, architecture, halo_size
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: inactive_node
using Oceananigans.Utils
using Adapt
using KernelAbstractions: @kernel, @index

## The equations are solved in an iterative form following the EVP rheology of
## Kimmritz et al (2016) (https://www.sciencedirect.com/science/article/pii/S1463500317300690)
#
# Where:
# σᵢⱼ(u) = 2η ϵ̇ᵢⱼ + [(ζ - η) * (ϵ̇₁₁ + ϵ̇₂₂) - P / 2] δᵢⱼ
#
struct ElastoViscoPlasticRheology{FT}
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    yield_curve_eccentricity :: FT # elliptic yield curve eccentricity
    minimum_plastic_stress :: FT # minimum plastic parameter (transitions to viscous behaviour)
    min_relaxation_parameter :: FT # minimum number of substeps expressed as the dynamic coefficient
    max_relaxation_parameter :: FT # maximum number of substeps expressed as the dynamic coefficient
end

"""
    ElastoViscoPlasticRheology(FT::DataType = Float64; 
                               ice_compressive_strength = 27500, 
                               ice_compaction_hardening = 20, 
                               yield_curve_eccentricity = 2, 
                               minimum_plastic_stress = 2e-9,
                               min_relaxation_parameter = 50,
                               max_relaxation_parameter = 300)

Constructs an `ElastoViscoPlasticRheology` object representing a "modified" elasto-visco-plastic
rheology for slab sea ice dynamics that follows the implementation of Kimmritz et al (2016).
The stress tensor is computed following the constitutive relation:
```math
σᵢⱼ = 2η ϵ̇ᵢⱼ + [(ζ - η) * (ϵ̇₁₁ + ϵ̇₂₂) - P / 2] δᵢⱼ
```
where ``ϵ̇ᵢⱼ`` are the strain rates, ``η`` is the shear viscosity, ``ζ`` is the bulk viscosity,
and ``P`` is the ice strength (acting as the isotropic part of the stress tensor)
parameterized as ``P★ h exp( - C ⋅ ( 1 - ℵ ))`` where ``P★`` is the `ice_compressive_strength`, 
``C`` is the `ice_compaction_hardening`, ``h`` is the ice thickness, and ``ℵ`` is the ice concentration.

The stresses are substepped using a dynamic substepping coefficient ``α`` that is
spatially varying and computed dynamically as in Kimmritz et al (2016)
In particular: α = sqrt(γ²) 
where γ² = ζ * π² * (Δt / mᵢ) / Az is a stability parameter with ``Az`` is the area of the grid cell, 
``mᵢ`` the ice mass, and ``Δt`` the time step.

The stresses are substepped with:
```math
σᵢⱼᵖ⁺¹ = σᵢⱼᵖ + (σᵢⱼᵖ⁺¹ - σᵢⱼᵖ) / α
```

This formulation allows fast convergence in regions where α is small. Regions where
α is large correspond to regions where the ice is more solid and the convergence is slower.
α can be thougth of as a ``pseudo substep number'' or a ``relaxation parameter''. If we are using 
a subcycling solver, if `α` ≪ number of substeps, the convergence will be faster.

Arguments
=========
    
- `grid`: the `SlabSeaIceModel` grid

Keyword Arguments
=================
    
- `ice_compressive_strength`: parameter expressing compressive strength (in Nm²). Default `27500`.
- `ice_compaction_hardening`: exponent coefficient for compaction hardening. Default `20`.
- `yield_curve_eccentricity`: eccentricity of the elliptic yield curve. Default `2`.
- `Δ_min`: Minimum value for the visco-plastic parameter. Limits the maximum viscosity of the ice, 
           transitioning the ice from a plastic to a viscous behaviour. Default value is `1e-10`.
- `min_relaxation_parameter`: Minimum value for the relaxation parameter `α`. Default value is `30`.
- `max_relaxation_parameter`: Maximum value for the relaxation parameter `α`. Default value is `500`.
"""
function ElastoViscoPlasticRheology(FT::DataType = Float64; 
                                    ice_compressive_strength = 27500, 
                                    ice_compaction_hardening = 20, 
                                    yield_curve_eccentricity = 2, 
                                    minimum_plastic_stress = 2e-9,
                                    min_relaxation_parameter = 50,
                                    max_relaxation_parameter = 300)

    return ElastoViscoPlasticRheology(convert(FT, ice_compressive_strength), 
                                      convert(FT, ice_compaction_hardening), 
                                      convert(FT, yield_curve_eccentricity),
                                      convert(FT, minimum_plastic_stress),
                                      convert(FT, min_relaxation_parameter),
                                      convert(FT, max_relaxation_parameter))
end

function rheology_auxiliary_fields(r::ElastoViscoPlasticRheology, grid)
    
    # TODO: What about boundary conditions?
    σ₁₁ = Field{Center, Center, Nothing}(grid)
    σ₂₂ = Field{Center, Center, Nothing}(grid)
    σ₁₂ = Field{Center, Center, Nothing}(grid)

    P  = Field{Center, Center, Nothing}(grid)
    ζ  = Field{Center, Center, Nothing}(grid)
    Δ  = Field{Center, Center, Nothing}(grid)

    α  = Field{Face, Face, Nothing}(grid) # Dynamic substeps a la Kimmritz et al (2016)
    uⁿ = Field{Face, Face, Nothing}(grid)
    vⁿ = Field{Face, Face, Nothing}(grid)

    # An initial (safe) educated guess
    fill!(α, r.max_relaxation_parameter)

    return (; σ₁₁, σ₂₂, σ₁₂, ζ, Δ, α, uⁿ, vⁿ, P)
end

# Extend the `adapt_structure` function for the ElastoViscoPlasticRheology
Adapt.adapt_structure(to, r::ElastoViscoPlasticRheology) = 
    ElastoViscoPlasticRheology(Adapt.adapt(to, r.ice_compressive_strength),
                               Adapt.adapt(to, r.ice_compaction_hardening),
                               Adapt.adapt(to, r.yield_curve_eccentricity),
                               Adapt.adapt(to, r.minimum_plastic_stress),
                               Adapt.adapt(to, r.min_relaxation_parameter),
                               Adapt.adapt(to, r.max_relaxation_parameter))

"""
    initialize_rheology!(model, rheology::ElastoViscoPlasticRheology)

Initialize the elasto-visco-plastic rheology.
In this step we calculate the ice strength given the ice mass (thickness and concentration).
"""
function initialize_rheology!(model, rheology::ElastoViscoPlasticRheology)
    h = model.ice_thickness
    ℵ = model.ice_concentration

    P★ = rheology.ice_compressive_strength
    C  = rheology.ice_compaction_hardening
    
    u, v   = model.velocities
    fields = model.dynamics.auxiliary_fields

    # compute on the whole grid including halos
    parameters = KernelParameters(size(fields.P.data)[1:2], fields.P.data.offsets[1:2])
    launch!(architecture(model.grid), model.grid, parameters, _initialize_evp_rhology!, fields, model.grid, P★, C, h, ℵ, u, v)
    
    return nothing
end

@kernel function _initialize_evp_rhology!(fields, grid, P★, C, h, ℵ, u, v)
    i, j = @index(Global, NTuple)    
    @inbounds fields.P[i, j, 1]  = ice_strength(i, j, 1, grid, P★, C, h, ℵ)
    @inbounds fields.uⁿ[i, j, 1] = u[i, j, 1]
    @inbounds fields.vⁿ[i, j, 1] = v[i, j, 1]
end

# The parameterization for an `ElastoViscoPlasticRheology`
@inline ice_strength(i, j, k, grid, P★, C, h, ℵ) = @inbounds P★ * h[i, j, k] * exp(- C * (1 - ℵ[i, j, k])) 

# Specific compute stresses for the EVP rheology
function compute_stresses!(model, dynamics, rheology::ElastoViscoPlasticRheology, Δt, Ns) 

    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ρᵢ = model.ice_density
    ℵ  = model.ice_concentration

    fields = dynamics.auxiliary_fields
    u, v = model.velocities

    Nx, Ny, _ = size(grid)
    Hx, Hy, _ = halo_size(grid)

    parameters = KernelParameters(-Hx+2:Nx+Hx-1, -Hy+2:Ny+Hy-1)

    launch!(arch, grid, parameters, _compute_evp_viscosities!, fields, grid, rheology, u, v, h, ℵ, ρᵢ, Δt)

    return nothing
end

@inline strain_rate_xx(i, j, k, grid, u, v) =  ℑyᵃᶜᵃ(i, j, k, grid, δxᶜᵃᵃ, Δy_qᶠᶠᶜ, u) / Azᶜᶜᶜ(i, j, k, grid)
@inline strain_rate_yy(i, j, k, grid, u, v) =  ℑxᶜᵃᵃ(i, j, k, grid, δyᵃᶜᵃ, Δx_qᶠᶠᶜ, v) / Azᶜᶜᶜ(i, j, k, grid)
@inline strain_rate_xy(i, j, k, grid, u, v) = (ℑyᵃᶜᵃ(i, j, k, grid, δxᶜᵃᵃ, Δy_qᶠᶠᶜ, v) + ℑxᶜᵃᵃ(i, j, k, grid, δyᵃᶜᵃ, Δx_qᶠᶠᶜ, u)) / Azᶜᶜᶜ(i, j, k, grid) / 2

@kernel function _compute_evp_viscosities!(fields, grid, rheology, u, v, h, ℵ, ρᵢ, Δt)
    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.minimum_plastic_stress

    # Extract auxiliary fields 
    P = fields.P

    # Strain rates
    ϵ̇₁₁ = strain_rate_xx(i, j, kᴺ, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, kᴺ, grid, u, v) 
    ϵ̇₁₂ = strain_rate_xy(i, j, kᴺ, grid, u, v)

    # Ice divergence 
    δ = ϵ̇₁₁ + ϵ̇₂₂

    # Ice shear (at Centers)
    s = sqrt((ϵ̇₁₁ - ϵ̇₂₂)^2 + 4ϵ̇₁₂^2)

    # Visco - Plastic parameter 
    # if Δ is very small we assume a linear viscous response
    # adding a minimum Δ_min (at Centers)
    Δ = max(sqrt(δ^2 + s^2 * e⁻²), Δm)
    P = @inbounds P[i, j, kᴺ]
    ζ = P / 2Δ

    @inbounds fields.ζ[i, j, kᴺ] = ζ
    @inbounds fields.Δ[i, j, kᴺ] = Δ

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.minimum_plastic_stress
    α⁺  = rheology.max_relaxation_parameter
    α⁻  = rheology.min_relaxation_parameter

    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂
    α   = fields.α
    
    # replacement pressure?
    Pᵣ = P * Δ / (Δ + Δm)
    η  = ζ * e⁻²

    # σ(uᵖ): the tangential stress depends only shear viscosity 
    # while the compressive stresses depend on the bulk viscosity and the ice strength
    σ₁₁ᵖ⁺¹ = 2 * η * ϵ̇₁₁ + ((ζ - η) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2) 
    σ₂₂ᵖ⁺¹ = 2 * η * ϵ̇₂₂ + ((ζ - η) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2)
    σ₁₂ᵖ⁺¹ = 2 * η * ϵ̇₁₂

    mᵢ = ice_mass(i, j, kᴺ, grid, h, ℵ, ρᵢ) 
    
    # Update coefficients for substepping using dynamic substepping
    # with spatially varying coefficients as in Kimmritz et al (2016)
    γ = ζ * π^2 * Δt / mᵢ / Azᶜᶜᶜ(i, j, kᴺ, grid)
    α = clamp(sqrt(γ), α⁻, α⁺)
    α = ifelse(isnan(α), α⁺, α)

    @inbounds begin
        # Compute the new stresses and store the value of the 
        # dynamic substepping coefficient α
        σ₁₁★ = (σ₁₁ᵖ⁺¹ - σ₁₁[i, j, kᴺ]) / α
        σ₂₂★ = (σ₂₂ᵖ⁺¹ - σ₂₂[i, j, kᴺ]) / α
        σ₁₂★ = (σ₁₂ᵖ⁺¹ - σ₁₂[i, j, kᴺ]) / α

        σ₁₁[i, j, kᴺ] += ifelse(mᵢ > 0, σ₁₁★, zero(grid))
        σ₂₂[i, j, kᴺ] += ifelse(mᵢ > 0, σ₂₂★, zero(grid))
        σ₁₂[i, j, kᴺ] += ifelse(mᵢ > 0, σ₁₂★, zero(grid))
        
        fields.α[i, j, kᴺ] = α 
    end
end

#####
##### Internal stress divergence for the EVP model
#####

# Here we extend all the functions that a rheology model needs to support:
@inline ice_stress_ux(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = ℑyᵃᶠᵃ(i, j, k, grid, fields.σ₁₁)
@inline ice_stress_uy(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = ℑxᶠᵃᵃ(i, j, k, grid, fields.σ₁₂)
@inline ice_stress_vx(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = ℑyᵃᶠᵃ(i, j, k, grid, fields.σ₁₂)
@inline ice_stress_vy(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = ℑxᶠᵃᵃ(i, j, k, grid, fields.σ₂₂)

# To help convergence to the right velocities
@inline compute_substep_Δtᶠᶠᶜ(i, j, grid, Δt, ::ElastoViscoPlasticRheology, substeps, fields) = @inbounds Δt / ℑxyᶠᶠᵃ(i, j, 1, grid, fields.α)

#####
##### Numerical forcing to help convergence
#####

@inline function sum_of_forcing_u(i, j, k, grid, ::ElastoViscoPlasticRheology, u_forcing, fields, Δt) 
    user_forcing = u_forcing(i, j, k, grid, fields)
    rheology_forcing = @inbounds (fields.uⁿ[i, j, k] - fields.u[i, j, k]) / Δt / ℑxyᶠᶠᵃ(i, j, k, grid, fields.α)
    return user_forcing + rheology_forcing
end

@inline function sum_of_forcing_v(i, j, k, grid, ::ElastoViscoPlasticRheology, v_forcing, fields, Δt) 
    user_forcing = v_forcing(i, j, k, grid, fields)
    rheology_forcing = @inbounds (fields.vⁿ[i, j, k] - fields.v[i, j, k]) / Δt / ℑxyᶠᶠᵃ(i, j, k, grid, fields.α)
    return user_forcing + rheology_forcing
end