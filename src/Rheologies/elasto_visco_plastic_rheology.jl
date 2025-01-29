using Oceananigans.TimeSteppers: store_field_tendencies!
using Oceananigans.Operators
using Oceananigans.Grids: AbstractGrid, architecture
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
    min_substeps :: FT # minimum number of substeps expressed as the dynamic coefficient
    max_substeps :: FT # maximum number of substeps expressed as the dynamic coefficient
end

"""
    ElastoViscoPlasticRheology(grid::AbstractGrid; 
                                 ice_compressive_strength = 27500, 
                                 ice_compaction_hardening = 20, 
                                 yield_curve_eccentricity = 2, 
                                 Δ_min = 1e-10,
                                 min_substeps = 30,
                                 max_substeps = 500)

Constructs an `ElastoViscoPlasticRheology` object representing a "modified" elasto-visco-plastic
rheology for slab sea ice dynamics that follows the implementation of Kimmritz et al (2016).
The `ice_strength` is parameterized as ``P★ h exp( - C ⋅ ( 1 - ℵ ))``. The number of substeps is
spatially varying and computed dynamically as in Kimmritz et al (2016):

In particular: α (stress substeps) == β (momentum substeps) = sqrt(γ) 
where γ = ζ * π² * (Δt / mᵢ) / Az is a stability parameter with ``Az`` the area of the grid cell, 
``ζ`` the bulk viscosity, ``mᵢ`` the ice mass, and ``Δt`` the time step.

This formulation allows fast convergence in regions where sqrt(γ) is small. Regions where
the coefficients are large correspond to regions where the ice is more solid and moves as a block
and the convergence is slower.

The number of substeps is then bounded by `min_substeps` and `max_substeps`.

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
- `min_substeps`: Minimum number of substeps expressed as the dynamic coefficient. Default value is `30`.
- `max_substeps`: Maximum number of substeps expressed as the dynamic coefficient. Default value is `500`.
"""
function ElastoViscoPlasticRheology(FT::DataType = Float64; 
                                    ice_compressive_strength = 27500, 
                                    ice_compaction_hardening = 20, 
                                    yield_curve_eccentricity = 2, 
                                    minimum_plastic_stress = 2e-9,
                                    min_substeps = 30,
                                    max_substeps = 1000)

    return ElastoViscoPlasticRheology(convert(FT, ice_compressive_strength), 
                                      convert(FT, ice_compaction_hardening), 
                                      convert(FT, yield_curve_eccentricity),
                                      convert(FT, minimum_plastic_stress),
                                      convert(FT, min_substeps),
                                      convert(FT, max_substeps))
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
    fill!(α, r.max_substeps)

    return (; σ₁₁, σ₂₂, σ₁₂, ζ, Δ, α, uⁿ, vⁿ, P)
end

# Extend the `adapt_structure` function for the ElastoViscoPlasticRheology
Adapt.adapt_structure(to, r::ElastoViscoPlasticRheology) = 
    ElastoViscoPlasticRheology(Adapt.adapt(to, r.ice_compressive_strength),
                               Adapt.adapt(to, r.ice_compaction_hardening),
                               Adapt.adapt(to, r.yield_curve_eccentricity),
                               Adapt.adapt(to, r.minimum_plastic_stress),
                               Adapt.adapt(to, r.min_substeps),
                               Adapt.adapt(to, r.max_substeps))

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
    fields = model.ice_dynamics.auxiliary_fields

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
function compute_stresses!(model, ice_dynamics, rheology::ElastoViscoPlasticRheology, Δt) 

    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ρᵢ = model.ice_density
    ℵ  = model.ice_concentration

    fields = ice_dynamics.auxiliary_fields
    u, v = model.velocities

    Nx, Ny, _ = size(grid)

    parameters = KernelParameters(-1:Nx+2, -1:Ny+2)

    launch!(arch, grid, parameters, _compute_evp_viscosities!, fields, grid, rheology, u, v)
    launch!(arch, grid, parameters, _compute_evp_stresses!, fields, grid, rheology, u, v, h, ℵ, ρᵢ, Δt)

    return nothing
end

const c = Center()
const f = Face()

# Hardcode No-slip boundary conditions on immersed boundaries?
# TODO: find a way not to hard-code. We need to pass the immersed_boundary_conditions of
# the velocities to the kernels
@inline function ∂xᴮᶜᶜᶜ(i, j, k, grid, u)
    i1 = inactive_node(i,   j, k, grid, f, c, c)
    i2 = inactive_node(i+1, j, k, grid, f, c, c) 
    Δx = Δxᶜᶜᶜ(i, j, k, grid)

    u1 = @inbounds u[i,   j, k]
    u2 = @inbounds u[i+1, j, k]

    return ifelse(i1, 2u2 / Δx, ifelse(i2, - 2u1 / Δx, (u2 - u1) / Δx))
end

@inline function ∂yᴮᶜᶜᶜ(i, j, k, grid, v)
    j1 = inactive_node(i, j,   k, grid, c, f, c)
    j2 = inactive_node(i, j+1, k, grid, c, f, c) 
    Δy = Δyᶜᶜᶜ(i, j, k, grid)

    v1 = @inbounds v[i, j,   k]
    v2 = @inbounds v[i, j+1, k]

    return ifelse(j1, 2v2 / Δy, ifelse(j2, 2v1 / Δy, (v2 - v1) / Δy))
end

@inline function ∂xᴮᶠᶠᶜ(i, j, k, grid, v)
    i1 = inactive_node(i-1, j, k, grid, c, f, c)
    i2 = inactive_node(i,   j, k, grid, c, f, c) 
    Δx = Δxᶠᶠᶜ(i, j, k, grid)

    v1 = @inbounds v[i-1, j, k]
    v2 = @inbounds v[i,   j, k]

    return ifelse(i1, 2v2 / Δx, ifelse(i2, - 2v1 / Δx, (v2 - v1) / Δx))
end

@inline function ∂yᴮᶠᶠᶜ(i, j, k, grid, u)
    j1 = inactive_node(i, j-1, k, grid, f, c, c)
    j2 = inactive_node(i, j,   k, grid, f, c, c) 
    Δy = Δyᶜᶜᶜ(i, j, k, grid)

    u1 = @inbounds u[i, j-1, k]
    u2 = @inbounds u[i, j,   k]

    return ifelse(j1, 2u2 / Δy, ifelse(j2, 2u1 / Δy, (u2 - u1) / Δy))
end

@kernel function _compute_evp_viscosities!(fields, grid, rheology, u, v)
    i, j = @index(Global, NTuple)

    P = fields.P

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.minimum_plastic_stress

    # Extract auxiliary fields 
    P = fields.P

    # Strain rates
    ϵ̇₁₁ =  ∂xᴮᶜᶜᶜ(i, j, 1, grid, u)
    ϵ̇₂₂ =  ∂yᴮᶜᶜᶜ(i, j, 1, grid, v)

    # Center - Center variables:
    ϵ̇₁₂ᶜᶜᶜ = (ℑxyᶜᶜᵃ(i, j, 1, grid, ∂xᴮᶠᶠᶜ, v) + 
              ℑxyᶜᶜᵃ(i, j, 1, grid, ∂yᴮᶠᶠᶜ, u)) / 2

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

# Compute the visco-plastic stresses for a slab sea ice model.
# The function updates the internal stress variables `σ₁₁`, `σ₂₂`, and `σ₁₂` in the `rheology` object
# following the mEVP formulation of Kimmritz et al (2016).
# This is the `meat` of the formulation.
@kernel function _compute_evp_stresses!(fields, grid, rheology, u, v, h, ℵ, ρᵢ, Δt)
    i, j = @index(Global, NTuple)

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.minimum_plastic_stress
    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂
    α   = fields.α

    # Strain rates
    ϵ̇₁₁ =  ∂xᴮᶜᶜᶜ(i, j, 1, grid, u)
    ϵ̇₁₂ = (∂xᴮᶠᶠᶜ(i, j, 1, grid, v) + ∂yᴮᶠᶠᶜ(i, j, 1, grid, u)) / 2
    ϵ̇₂₂ =  ∂yᴮᶜᶜᶜ(i, j, 1, grid, v)

    Pᶜᶜᶜ = @inbounds fields.P[i, j, 1]
    ζᶜᶜᶜ = @inbounds fields.ζ[i, j, 1]
    Δᶜᶜᶜ = @inbounds fields.Δ[i, j, 1]
    ζᶠᶠᶜ = ℑxyᶠᶠᵃ(i, j, 1, grid, fields.ζ)

    # replacement pressure?
    Pᵣ = Pᶜᶜᶜ * Δᶜᶜᶜ / (Δᶜᶜᶜ + Δm)

    ηᶜᶜᶜ = ζᶜᶜᶜ * e⁻²
    ηᶠᶠᶜ = ζᶠᶠᶜ * e⁻²

    # σ(uᵖ): the tangential stress depends only shear viscosity 
    # while the compressive stresses depend on the bulk viscosity and the ice strength
    σ₁₁ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₁₁ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2) 
    σ₂₂ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₂₂ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2)
    σ₁₂ᵖ⁺¹ = 2 * ηᶠᶠᶜ * ϵ̇₁₂

    mᵢᶜᶜᶜ = ice_mass(i, j, 1, grid, h, ℵ, ρᵢ) 
    mᵢᶠᶠᶜ = ℑxyᶠᶠᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ) 

    # Update coefficients for substepping if we are using dynamic substepping
    # with spatially varying coefficients such as in Kimmritz et al (2016)
    γ²ᶜᶜᶜ = ζᶜᶜᶜ * π^2 * Δt / mᵢᶜᶜᶜ / Azᶜᶜᶜ(i, j, 1, grid)
    γᶜᶜᶜ  = clamp(sqrt(γ²ᶜᶜᶜ), rheology.min_substeps, rheology.max_substeps)
    γᶜᶜᶜ  = ifelse(isnan(γᶜᶜᶜ), rheology.max_substeps, γᶜᶜᶜ) # In case both ζᶜᶜᶜ and mᵢᶜᶜᶜ are zero

    γ²ᶠᶠᶜ = ζᶠᶠᶜ * π^2 * Δt / mᵢᶠᶠᶜ / Azᶠᶠᶜ(i, j, 1, grid)
    γᶠᶠᶜ  = clamp(sqrt(γ²ᶠᶠᶜ), rheology.min_substeps, rheology.max_substeps)
    γᶠᶠᶜ  = ifelse(isnan(γᶠᶠᶜ), rheology.max_substeps, γᶠᶠᶜ) # In case both ζᶠᶠᶜ and mᵢᶠᶠᶜ are zero

    # Compute the new stresses and store the value of the dynamic substepping coefficient α
    @inbounds σ₁₁[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, (σ₁₁ᵖ⁺¹ - σ₁₁[i, j, 1]) / γᶜᶜᶜ, zero(grid))
    @inbounds σ₂₂[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, (σ₂₂ᵖ⁺¹ - σ₂₂[i, j, 1]) / γᶜᶜᶜ, zero(grid))
    @inbounds σ₁₂[i, j, 1] += ifelse(mᵢᶠᶠᶜ > 0, (σ₁₂ᵖ⁺¹ - σ₁₂[i, j, 1]) / γᶠᶠᶜ, zero(grid))
    @inbounds   α[i, j, 1]  = γᶜᶜᶜ
    
    # Mask inactive nodes
    @inbounds σ₁₁[i, j, 1] = ifelse(inactive_node(i, j, 1, grid, Center(), Center(), Center()), zero(grid), σ₁₁[i, j, 1])
    @inbounds σ₂₂[i, j, 1] = ifelse(inactive_node(i, j, 1, grid, Center(), Center(), Center()), zero(grid), σ₂₂[i, j, 1])
    @inbounds σ₁₂[i, j, 1] = ifelse(inactive_node(i, j, 1, grid, Face(),   Face(),   Center()), zero(grid), σ₁₂[i, j, 1]) 
end

#####
##### Internal stress divergence for the EVP model
#####

# Here we extend all the functions that a rheology model needs to support:
@inline function ∂ⱼ_σ₁ⱼ(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields, Δt) 
    ∂xσ₁₁ = ∂xᶠᶜᶜ(i, j, k, grid, fields.σ₁₁)
    ∂yσ₁₂ = ∂yᶠᶜᶜ(i, j, k, grid, fields.σ₁₂)

    return ∂xσ₁₁ + ∂yσ₁₂
end

@inline function ∂ⱼ_σ₂ⱼ(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields, Δt) 
    ∂xσ₁₂ = ∂xᶜᶠᶜ(i, j, k, grid, fields.σ₁₂)
    ∂yσ₂₂ = ∂yᶜᶠᶜ(i, j, k, grid, fields.σ₂₂)

    return ∂xσ₁₂ + ∂yσ₂₂
end

# To help convergence to the right velocities
@inline compute_time_stepᶠᶜᶜ(i, j, grid, Δt, ::ElastoViscoPlasticRheology, substeps, fields) = Δt / ℑxᶠᵃᵃ(i, j, 1, grid, fields.α)
@inline compute_time_stepᶜᶠᶜ(i, j, grid, Δt, ::ElastoViscoPlasticRheology, substeps, fields) = Δt / ℑyᵃᶠᵃ(i, j, 1, grid, fields.α)

#####
##### Numerical forcing to help convergence
#####

@inline function sum_of_forcing_x(i, j, k, grid, ::ElastoViscoPlasticRheology, u_forcing, fields, Δt) 
    user_forcing = u_forcing(i, j, k, grid, fields)
    rheology_forcing = @inbounds (fields.uⁿ[i, j, k] - fields.u[i, j, k]) / Δt / ℑxᶠᵃᵃ(i, j, k, grid, fields.α)
    return user_forcing + rheology_forcing
end

@inline function sum_of_forcing_y(i, j, k, grid, ::ElastoViscoPlasticRheology, v_forcing, fields, Δt) 
    user_forcing = v_forcing(i, j, k, grid, fields)
    rheology_forcing = @inbounds (fields.vⁿ[i, j, k] - fields.v[i, j, k]) / Δt / ℑyᵃᶠᵃ(i, j, k, grid, fields.α)
    return user_forcing + rheology_forcing
end