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
struct ElastoViscoPlasticRheology{FT, IP}
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    yield_curve_eccentricity :: FT # elliptic yield curve eccentricity
    minimum_plastic_stress :: FT # minimum plastic parameter (transitions to viscous behaviour)
    min_relaxation_parameter :: FT # minimum number of substeps expressed as the dynamic coefficient
    max_relaxation_parameter :: FT # maximum number of substeps expressed as the dynamic coefficient
    relaxation_strength :: FT # strength of the relaxation parameter
    pressure_formulation :: IP # formulation of ice pressure
    ElastoViscoPlasticRheology(P::FT, C::FT, e::FT, Δ_min::FT, α⁻::FT, α⁺::FT, c::FT, ip::IP) where {FT, IP}  = 
        new{FT, IP}(P, C, e, Δ_min, α⁻, α⁺, c, ip)
end

struct ReplacementPressure end
struct IceStrength end

"""
    ElastoViscoPlasticRheology(FT::DataType = Float64; 
                               ice_compressive_strength = 27500, 
                               ice_compaction_hardening = 20, 
                               yield_curve_eccentricity = 2, 
                               minimum_plastic_stress = 2e-9,
                               min_relaxation_parameter = 50,
                               max_relaxation_parameter = 300,
                               relaxation_strength = π^2,
                               pressure_formulation = ReplacementPressure())

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
where γ² = ζ * cα * (Δt / mᵢ) / Az is a stability parameter with ``Az`` is the area of the grid cell, 
``mᵢ`` the ice mass, ``Δt`` the time step, and ``cα`` a numerical stability parameter which controls the 
stregth of ``γ²``.

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
    
- `ice_compressive_strength`: parameter expressing compressive strength (in Nm²). Default: `27500`.
- `ice_compaction_hardening`: exponent coefficient for compaction hardening. Default: `20`.
- `yield_curve_eccentricity`: eccentricity of the elliptic yield curve. Default: `2`.
- `Δ_min`: Minimum value for the visco-plastic parameter. Limits the maximum viscosity of the ice, 
           transitioning the ice from a plastic to a viscous behaviour. Default: `1e-10`.
- `min_relaxation_parameter`: Minimum value for the relaxation parameter `α`. Default: `30`.
- `max_relaxation_parameter`: Maximum value for the relaxation parameter `α`. Default: `500`.
- `relaxation_strength`: parameter controlling the strength of the relaxation parameter. The maximum value is `π²`, see Kimmritz et al (2016). Default: `π² / 2`.
- `pressure_formulation`: can use `ReplacementPressure` or `IceStrength`. The replacement pressure formulation avoids ice motion in the absence of forcing. Default: `ReplacementPressure`.
"""
function ElastoViscoPlasticRheology(FT::DataType = Float64; 
                                    ice_compressive_strength = 27500, 
                                    ice_compaction_hardening = 20, 
                                    yield_curve_eccentricity = 2, 
                                    minimum_plastic_stress = 2e-9,
                                    min_relaxation_parameter = 50,
                                    max_relaxation_parameter = 300,
                                    relaxation_strength = π^2,
                                    pressure_formulation = ReplacementPressure())

    return ElastoViscoPlasticRheology(convert(FT, ice_compressive_strength), 
                                      convert(FT, ice_compaction_hardening), 
                                      convert(FT, yield_curve_eccentricity),
                                      convert(FT, minimum_plastic_stress),
                                      convert(FT, min_relaxation_parameter),
                                      convert(FT, max_relaxation_parameter),
                                      convert(FT, relaxation_strength),
                                      pressure_formulation)
end

function required_auxiliary_fields(r::ElastoViscoPlasticRheology, grid)
    
    # TODO: What about boundary conditions?
    σ₁₁ = Field{Center, Center, Nothing}(grid)
    σ₂₂ = Field{Center, Center, Nothing}(grid)
    σ₁₂ = Field{Face, Face, Nothing}(grid)

    uⁿ = Field{Face,   Center, Nothing}(grid)
    vⁿ = Field{Center, Face,   Nothing}(grid)
    P  = Field{Center, Center, Nothing}(grid)
    α  = Field{Center, Center, Nothing}(grid) # Dynamic substeps a la Kimmritz et al (2016)
    ζ  = Field{Center, Center, Nothing}(grid)
    Δ  = Field{Center, Center, Nothing}(grid)

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
                               Adapt.adapt(to, r.max_relaxation_parameter),
                               Adapt.adapt(to, r.relaxation_strength),
                               Adapt.adapt(to, r.pressure_formulation))

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
function compute_stresses!(model, dynamics, rheology::ElastoViscoPlasticRheology, Δt) 

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

    launch!(arch, grid, parameters, _compute_evp_viscosities!, fields, grid, rheology, u, v)
    launch!(arch, grid, parameters, _compute_evp_stresses!, fields, grid, rheology, u, v, h, ℵ, ρᵢ, Δt)

    return nothing
end

@inline strain_rate_xx(i, j, k, grid, u, v) =  δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, u) / Azᶜᶜᶜ(i, j, k, grid)
@inline strain_rate_yy(i, j, k, grid, u, v) =  δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, v) / Azᶜᶜᶜ(i, j, k, grid)
@inline strain_rate_xy(i, j, k, grid, u, v) = (δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶠᶜ, v) + δyᵃᶠᵃ(i, j, k, grid, Δx_qᶠᶜᶜ, u)) / Azᶠᶠᶜ(i, j, k, grid) / 2

@kernel function _compute_evp_viscosities!(fields, grid, rheology, u, v)
    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.minimum_plastic_stress

    # Extract auxiliary fields 
    P = fields.P

    # Strain rates
    ϵ̇₁₁ = strain_rate_xx(i, j, kᴺ, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, kᴺ, grid, u, v) 

    # Center - Center variables:
    ϵ̇₁₂ᶜᶜᶜ = ℑxyᶜᶜᵃ(i, j, kᴺ, grid, strain_rate_xy, u, v)

    # Ice divergence 
    δ = ϵ̇₁₁ + ϵ̇₂₂

    # Ice shear (at Centers)
    s = sqrt((ϵ̇₁₁ - ϵ̇₂₂)^2 + 4ϵ̇₁₂ᶜᶜᶜ^2)

    # Visco - Plastic parameter 
    # if Δ is very small we assume a linear viscous response
    # adding a minimum Δ_min (at Centers)
    Δᶜᶜᶜ = max(sqrt(δ^2 + s^2 * e⁻²), Δm)
    Pᶜᶜᶜ = @inbounds P[i, j, 1]

    @inbounds fields.ζ[i, j, 1] = Pᶜᶜᶜ / 2Δᶜᶜᶜ
    @inbounds fields.Δ[i, j, 1] = Δᶜᶜᶜ
end

@inline ice_pressure(i, j, k, grid, ::IceStrength, r, fields) = @inbounds fields.P[i, j, k]

@inline function ice_pressure(i, j, k, grid, ::ReplacementPressure, r, fields)
    Pᶜᶜᶜ = @inbounds fields.P[i, j, k]
    Δᶜᶜᶜ = @inbounds fields.Δ[i, j, k]
    Δm   = r.minimum_plastic_stress
    return Pᶜᶜᶜ * Δᶜᶜᶜ / (Δᶜᶜᶜ + Δm)
end

# Compute the visco-plastic stresses for a slab sea ice model.
# The function updates the internal stress variables `σ₁₁`, `σ₂₂`, and `σ₁₂` in the `rheology` object
# following the αEVP formulation of Kimmritz et al (2016).
@kernel function _compute_evp_stresses!(fields, grid, rheology, u, v, h, ℵ, ρᵢ, Δt)
    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    α⁺  = rheology.max_relaxation_parameter
    α⁻  = rheology.min_relaxation_parameter
    cα  = rheology.relaxation_strength
    ip  = rheology.pressure_formulation

    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂
    α   = fields.α
    
    # Strain rates
    ϵ̇₁₁ = strain_rate_xx(i, j, kᴺ, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, kᴺ, grid, u, v) 
    ϵ̇₁₂ = strain_rate_xy(i, j, kᴺ, grid, u, v)

    ζᶜᶜᶜ = @inbounds fields.ζ[i, j, 1]
    ζᶠᶠᶜ = ℑxyᶠᶠᵃ(i, j, 1, grid, fields.ζ)

    # replacement pressure?
    Pᵣ = ice_pressure(i, j, 1, grid, ip, rheology, fields)

    ηᶜᶜᶜ = ζᶜᶜᶜ * e⁻²
    ηᶠᶠᶜ = ζᶠᶠᶜ * e⁻²

    # σ(uᵖ): the tangential stress depends only shear viscosity 
    # while the compressive stresses depend on the bulk viscosity and the ice strength
    σ₁₁ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₁₁ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2) 
    σ₂₂ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₂₂ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2)
    σ₁₂ᵖ⁺¹ = 2 * ηᶠᶠᶜ * ϵ̇₁₂

    mᵢᶜᶜᶜ = ice_mass(i, j, 1, grid, h, ℵ, ρᵢ) 
    mᵢᶠᶠᶜ = ℑxyᶠᶠᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ) 

    # Update coefficients for substepping using dynamic substepping
    # with spatially varying coefficients as in Kimmritz et al (2016)
    γ²ᶜᶜᶜ = ζᶜᶜᶜ * cα * Δt / mᵢᶜᶜᶜ / Azᶜᶜᶜ(i, j, 1, grid)
    γ²ᶜᶜᶜ = ifelse(isnan(γ²ᶜᶜᶜ), α⁺^2, γ²ᶜᶜᶜ) # In case both ζᶜᶜᶜ and mᵢᶜᶜᶜ are zero
    γᶜᶜᶜ  = clamp(sqrt(γ²ᶜᶜᶜ), α⁻, α⁺)

    γ²ᶠᶠᶜ = ζᶠᶠᶜ * cα * Δt / mᵢᶠᶠᶜ / Azᶠᶠᶜ(i, j, 1, grid)
    γ²ᶠᶠᶜ = ifelse(isnan(γ²ᶠᶠᶜ), α⁺^2, γ²ᶠᶠᶜ) # In case both ζᶠᶠᶜ and mᵢᶠᶠᶜ are zero
    γᶠᶠᶜ  = clamp(sqrt(γ²ᶠᶠᶜ), α⁻, α⁺)

    @inbounds begin
        # Compute the new stresses and store the value of the 
        # dynamic substepping coefficient α
        σ₁₁★ = (σ₁₁ᵖ⁺¹ - σ₁₁[i, j, 1]) / γᶜᶜᶜ
        σ₂₂★ = (σ₂₂ᵖ⁺¹ - σ₂₂[i, j, 1]) / γᶜᶜᶜ
        σ₁₂★ = (σ₁₂ᵖ⁺¹ - σ₁₂[i, j, 1]) / γᶠᶠᶜ

        σ₁₁[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, σ₁₁★, zero(grid))
        σ₂₂[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, σ₂₂★, zero(grid))
        σ₁₂[i, j, 1] += ifelse(mᵢᶠᶠᶜ > 0, σ₁₂★, zero(grid))
          α[i, j, 1]  = γᶜᶜᶜ
    end
end

#####
##### Internal stress divergence for the EVP model
#####

# Here we extend all the functions that a rheology model needs to support:
@inline ice_stress_ux(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = @inbounds fields.σ₁₁[i, j, k]
@inline ice_stress_vx(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = @inbounds fields.σ₁₂[i, j, k]
@inline ice_stress_uy(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = @inbounds fields.σ₁₂[i, j, k]
@inline ice_stress_vy(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = @inbounds fields.σ₂₂[i, j, k]

# To help convergence to the right velocities
@inline compute_substep_Δtᶠᶜᶜ(i, j, grid, Δt, ::ElastoViscoPlasticRheology, substeps, fields) = Δt / ℑxᶠᵃᵃ(i, j, 1, grid, fields.α)
@inline compute_substep_Δtᶜᶠᶜ(i, j, grid, Δt, ::ElastoViscoPlasticRheology, substeps, fields) = Δt / ℑyᵃᶠᵃ(i, j, 1, grid, fields.α)

#####
##### Numerical forcing to help convergence
#####

@inline function sum_of_forcing_u(i, j, k, grid, ::ElastoViscoPlasticRheology, u_forcing, fields, Δt) 
    user_forcing = u_forcing(i, j, k, grid, fields)
    rheology_forcing = @inbounds (fields.uⁿ[i, j, k] - fields.u[i, j, k]) / Δt / ℑxᶠᵃᵃ(i, j, k, grid, fields.α)
    return user_forcing + rheology_forcing
end

@inline function sum_of_forcing_v(i, j, k, grid, ::ElastoViscoPlasticRheology, v_forcing, fields, Δt) 
    user_forcing = v_forcing(i, j, k, grid, fields)
    rheology_forcing = @inbounds (fields.vⁿ[i, j, k] - fields.v[i, j, k]) / Δt / ℑyᵃᶠᵃ(i, j, k, grid, fields.α)
    return user_forcing + rheology_forcing
end
