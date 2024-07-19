using Oceananigans.TimeSteppers: store_field_tendencies!
using Oceananigans.Operators
using Oceananigans.Grids: AbstractGrid

## The equations are solved in an iterative form following the EVP rheology of
## Kimmritz et al (2016) (https://www.sciencedirect.com/science/article/pii/S1463500317300690)
#
# Where:
# σᵢⱼ(u) = 2η ϵ̇ᵢⱼ + [(ζ - η) * (ϵ̇₁₁ + ϵ̇₂₂) - P / 2] δᵢⱼ
#
struct ExplicitViscoPlasticRheology{FT} <: AbstractRheology
    cross_nudging :: FT
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    yield_curve_eccentricity :: FT # elliptic yield curve eccentricity
    Δ_min :: FT # minimum plastic parameter (transitions to viscous behaviour)
end

"""
    ExplicitViscoPlasticRheology(grid::AbstractGrid; 
                                 cross_nudging = 2,
                                 ice_compressive_strength = 27500, 
                                 ice_compaction_hardening = 20, 
                                 yield_curve_eccentricity = 2, 
                                 Δ_min = 2e-9)

Constructs an `ExplicitViscoPlasticRheology` object representing a "modified" elasto-visco-plastic
rheology for slab sea ice dynamics that follows the implementation of kimmritz et al (2016).
The `ice_strength` is parameterized as ``P★ h exp( - C ⋅ ( 1 - ℵ ))`` 

Arguments
=========
    
- `FT`: the floating point precision

Keyword Arguments
=================
    
- `ice_compressive_strength`: parameter expressing compressive strength (in Nm²), default `27500`.
- `ice_compaction_hardening`: exponent coefficient for compaction hardening, default `20`.
- `yield_curve_eccentricity`: eccentricity of the elliptic yield curve, default `2`.
- `Δ_min`: Minimum value for the visco-plastic parameter. Limits the maximum viscosity of the ice, 
           transitioning the ice from a plastic to a viscous behaviour. Default value is `2e-9`.
"""
function ExplicitViscoPlasticRheology(FT::DataType = Float64; 
                                      cross_nudging = 2,
                                      ice_compressive_strength = 27500, 
                                      ice_compaction_hardening = 20, 
                                      yield_curve_eccentricity = 2, 
                                      Δ_min = 2e-9)

    return ExplicitViscoPlasticRheology(convert(FT, cross_nudging),
                                        convert(FT, ice_compressive_strength), 
                                        convert(FT, ice_compaction_hardening), 
                                        convert(FT, yield_curve_eccentricity),
                                        convert(FT, Δ_min))
end

function required_auxiliary_fields(grid, ::ExplicitViscoPlasticRheology, ::CGridDynamics)

    σ₁₁ = Field{Center, Center, Nothing}(grid)
    σ₂₂ = Field{Center, Center, Nothing}(grid)
    σ₁₂ = Field{Face, Face, Nothing}(grid)

    uⁿ = Field{Face, Center, Nothing}(grid)
    vⁿ = Field{Center, Face, Nothing}(grid)
    P  = Field{Center, Center, Nothing}(grid)

    return (; σ₁₁, σ₂₂, σ₁₂, uⁿ, vⁿ, P)
end

function required_auxiliary_fields(grid, ::ExplicitViscoPlasticRheology, ::EGridDynamics)

    û = YFaceField(grid) 
    v̂ = XFaceField(grid) 

    σ₁₁ = CenterField(grid)
    σ₂₂ = CenterField(grid)
    σ₁₂ = Field{Face, Face, Center}(grid)

    σ̂₁₁ = Field{Face, Face, Center}(grid)
    σ̂₂₂ = Field{Face, Face, Center}(grid)
    σ̂₁₂ = CenterField(grid)

    uⁿ = XFaceField(grid) 
    vⁿ = YFaceField(grid) 

    ûⁿ = YFaceField(grid) 
    v̂ⁿ = XFaceField(grid) 

    P = Field{Face, Face, Center}(grid)
    P̂ = Field{Face, Face, Center}(grid)

    return (; û, v̂,
              σ₁₁, σ₂₂, σ₁₂, uⁿ, vⁿ, P,
              σ̂₁₁, σ̂₂₂, σ̂₁₂, ûⁿ, v̂ⁿ, P̂)
end

# Extend the `adapt_structure` function for the ExplicitViscoPlasticRheology
Adapt.adapt_structure(to, r::ExplicitViscoPlasticRheology) = 
    ExplicitViscoPlasticRheology(Adapt.adapt(to, r.cross_nudging),
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

    P★ = rheology.ice_compressive_strength
    C  = rheology.ice_compaction_hardening
    
    u, v   = model.velocities
    fields = model.ice_dynamics.auxiliary_fields
    dgrid  = dynamics_grid(model.ice_dynamics) 
    grid   = model.grid

    # compute on the whole grid including halos
    parameters = KernelParameters(size(fields.P.data)[1:2], fields.P.data.offsets[1:2])
    launch!(architecture(grid), grid, parameters, _initialize_evp_rhology!, fields, dgrid, grid, P★, C, h, ℵ, u, v)
    
    return nothing
end

# Initialize the rheology on a C-grid
@kernel function _initialize_evp_rhology!(fields, ::CGridDynamics, grid, P★, C, h, ℵ, u, v)
    i, j = @index(Global, NTuple)    
    @inbounds fields.P[i, j, 1]  = ice_strength(i, j, 1, grid, P★, C, h, ℵ)
    @inbounds fields.uⁿ[i, j, 1] = u[i, j, 1]
    @inbounds fields.vⁿ[i, j, 1] = v[i, j, 1]
end

# Initialize the rheology on an E-grid
@kernel function _initialize_evp_rhology!(fields, ::EGridDynamics, grid, P★, C, h, ℵ, u, v)
    i, j = @index(Global, NTuple)    
    @inbounds begin
        fields.P[i, j, 1]  = ice_strength(i, j, 1, grid, P★, C, h, ℵ)
        fields.P̂[i, j, 1]  = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, ice_strength, P★, C, h, ℵ)
        fields.uⁿ[i, j, 1] = u[i, j, 1]
        fields.vⁿ[i, j, 1] = v[i, j, 1]
        fields.ûⁿ[i, j, 1] = fields.û[i, j, 1]
        fields.v̂ⁿ[i, j, 1] = fields.v̂[i, j, 1]
    end
end

# The parameterization for an `ExplicitViscoPlasticRheology`
@inline ice_strength(i, j, k, grid, P★, C, h, ℵ) = P★ * h[i, j, k] * exp(- C * (1 - ℵ[i, j, k])) 

# Specific compute stresses for the EVP rheology
function compute_stresses!(model, solver, rheology::ExplicitViscoPlasticRheology, Δt) 

    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ℵ  = model.ice_concentration
    ρᵢ = model.ice_density

    fields   = solver.auxiliary_fields
    substeps = solver.substeps
    stepping_coefficient = solver.substepping_coefficient

    u, v = model.velocities
    _compute_stress_kernel! = dynamics_grid(solver) isa CGridDynamics ? _compute_cgrid_evp_stresses! : _compute_egrid_evp_stresses!

    launch!(arch, grid, :xyz, _compute_stress_kernel!, fields, rheology, grid, u, v, h, ℵ, ρᵢ, Δt, substeps, stepping_coefficient)

    return nothing
end

# Compute the visco-plastic stresses for a slab sea ice model.
# The function updates the internal stress variables `σ₁₁`, `σ₂₂`, and `σ₁₂` in the `rheology` object
# following the mEVP formulation of Kimmritz et al (2016).
# This is the `meat` of the formulation.
@kernel function _compute_cgrid_evp_stresses!(fields, rheology, grid, u, v, h, ℵ, ρᵢ, Δt, substeps, stepping_coefficient)
    i, j = @index(Global, NTuple)

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.Δ_min

    # Extract auxiliary fields 
    P   = fields.P
    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂

    # Strain rates
    ϵ̇₁₁ =  ∂xᶜᶜᶜ_U(i, j, 1, grid, u)
    ϵ̇₁₂ = (∂xᶠᶠᶜ_c(i, j, 1, grid, v) + ∂yᶠᶠᶜ_c(i, j, 1, grid, u)) / 2
    ϵ̇₂₂ =  ∂yᶜᶜᶜ_V(i, j, 1, grid, v)

    # Center - Center variables:
    ϵ̇₁₂ᶜᶜᶜ = (ℑxyᴮᶜᶜᶜ(i, j, 1, grid, ∂xᶠᶠᶜ_c, v) + ℑxyᴮᶜᶜᶜ(i, j, 1, grid, ∂yᶠᶠᶜ_c, u)) / 2

    # Ice divergence 
    δ = ϵ̇₁₁ + ϵ̇₂₂

    # Ice shear (at Centers)
    s = sqrt((ϵ̇₁₁ - ϵ̇₂₂)^2 + 4ϵ̇₁₂ᶜᶜᶜ^2)

    # Visco - Plastic parameter 
    # if Δ is very small we assume a linear viscous response
    # adding a minimum Δ_min (at Centers)
    Δᶜᶜᶜ = sqrt(δ^2 + s^2 * e⁻²) + Δm

    # Face - Face variables
    ϵ̇₁₁ᶠᶠᶜ = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, ∂xᶜᶜᶜ_U, u)
    ϵ̇₂₂ᶠᶠᶜ = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, ∂yᶜᶜᶜ_V, v)

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

    mᵢᶜᶜᶜ = ice_mass(i, j, 1, grid, h, ℵ, ρᵢ) 
    mᵢᶠᶠᶜ = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ) 

    c = stepping_coefficient

    # Update coefficients for substepping if we are using dynamic substepping
    # with spatially varying coefficients such as in Kimmritz et al (2016)
    update_stepping_coefficients!(i, j, 1, grid, c, ζᶜᶜᶜ, mᵢᶜᶜᶜ, Δt)

    # Coefficient for substepping internal stress
    αᶜᶜᶜ = get_stepping_coefficients(i, j, 1, grid, substeps, c)
    αᶠᶠᶜ = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, get_stepping_coefficients, substeps, c)

    @inbounds σ₁₁[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, (σ₁₁ᵖ⁺¹ - σ₁₁[i, j, 1]) / αᶜᶜᶜ, zero(grid))
    @inbounds σ₂₂[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, (σ₂₂ᵖ⁺¹ - σ₂₂[i, j, 1]) / αᶜᶜᶜ, zero(grid))
    @inbounds σ₁₂[i, j, 1] += ifelse(mᵢᶠᶠᶜ > 0, (σ₁₂ᵖ⁺¹ - σ₁₂[i, j, 1]) / αᶠᶠᶜ, zero(grid))
end

# Compute the visco-plastic stresses for a slab sea ice model.
# The function updates the internal stress variables `σ₁₁`, `σ₂₂`, and `σ₁₂` in the `rheology` object
# following the mEVP formulation of Kimmritz et al (2016).
# This is the `meat` of the formulation.
@kernel function _compute_egrid_evp_stresses!(fields, rheology, grid, u, v, h, ℵ, ρᵢ, Δt, substeps, stepping_coefficient)
    i, j = @index(Global, NTuple)

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.Δ_min
    γc  = rheology.cross_nudging

    # Extract auxiliary fields 
    P   = fields.P
    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂

    P̂   = fields.P̂
    σ̂₁₁ = fields.σ̂₁₁
    σ̂₂₂ = fields.σ̂₂₂
    σ̂₁₂ = fields.σ̂₁₂
    û   = fields.û
    v̂   = fields.v̂
 
    # Strain rates
    ϵ̇₁₁ =  ∂xᶜᶜᶜ_U(i, j, 1, grid, u)
    ϵ̇₁₂ = (∂xᶠᶠᶜ_c(i, j, 1, grid, v) + ∂yᶠᶠᶜ_c(i, j, 1, grid, u)) / 2
    ϵ̇₂₂ =  ∂yᶜᶜᶜ_V(i, j, 1, grid, v)

    # Complementary strain rates
    ϵ̇₁₁ᶠᶠᶜ =  ∂xᶠᶠᶜ_c(i, j, 1, grid, û)
    ϵ̇₁₂ᶜᶜᶜ = (∂xᶜᶜᶜ_U(i, j, 1, grid, v̂) + ∂yᶜᶜᶜ_V(i, j, 1, grid, û)) / 2
    ϵ̇₂₂ᶠᶠᶜ =  ∂yᶠᶠᶜ_c(i, j, 1, grid, v̂)

    # Ice divergence 
    δ = ϵ̇₁₁ + ϵ̇₂₂

    # Ice shear (at Centers)
    s = sqrt((ϵ̇₁₁ - ϵ̇₂₂)^2 + 4ϵ̇₁₂ᶜᶜᶜ^2)

    # Visco - Plastic parameter 
    # if Δ is very small we assume a linear viscous response
    # adding a minimum Δ_min (at Centers)
    Δᶜᶜᶜ = sqrt(δ^2 + s^2 * e⁻²) + Δm

    # Face - Face variables
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
    Pᶠᶠᶜ = @inbounds P̂[i, j, 1]

    # ζ: Bulk viscosity (viscosity which responds to compression) 
    # η: Shear viscosity (viscosity which responds to shear)
    ζᶜᶜᶜ = Pᶜᶜᶜ / (2Δᶜᶜᶜ)
    ηᶜᶜᶜ = ζᶜᶜᶜ * e⁻²

    ζᶠᶠᶜ = Pᶠᶠᶜ / (2Δᶠᶠᶜ)
    ηᶠᶠᶜ = ζᶠᶠᶜ * e⁻²

    # Replacement pressure
    Pᵣᶜᶜᶜ = Pᶜᶜᶜ * Δᶜᶜᶜ / (Δᶜᶜᶜ + Δm)
    Pᵣᶠᶠᶜ = Pᶠᶠᶜ * Δᶠᶠᶜ / (Δᶠᶠᶜ + Δm)

    # σ(uᵖ): the tangential stress depends only shear viscosity 
    # while the compressive stresses depend on the bulk viscosity and the ice strength
    σ₁₁ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₁₁ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣᶜᶜᶜ / 2) 
    σ₂₂ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₂₂ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣᶜᶜᶜ / 2)
    σ₁₂ᵖ⁺¹ = 2 * ηᶠᶠᶜ * ϵ̇₁₂

    σ̂₁₁ᵖ⁺¹ = 2 * ηᶠᶠᶜ * ϵ̇₁₁ᶠᶠᶜ + ((ζᶠᶠᶜ - ηᶠᶠᶜ) * (ϵ̇₁₁ᶠᶠᶜ + ϵ̇₂₂ᶠᶠᶜ) - Pᵣᶠᶠᶜ / 2) 
    σ̂₂₂ᵖ⁺¹ = 2 * ηᶠᶠᶜ * ϵ̇₂₂ᶠᶠᶜ + ((ζᶠᶠᶜ - ηᶠᶠᶜ) * (ϵ̇₁₁ᶠᶠᶜ + ϵ̇₂₂ᶠᶠᶜ) - Pᵣᶠᶠᶜ / 2)
    σ̂₁₂ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₁₂ᶜᶜᶜ

    mᵢᶜᶜᶜ = ice_mass(i, j, 1, grid, h, ℵ, ρᵢ) 
    mᵢᶠᶠᶜ = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ) 

    c = stepping_coefficient

    # Update coefficients for substepping if we are using dynamic substepping
    # with spatially varying coefficients such as in Kimmritz et al (2016)
    update_stepping_coefficients!(i, j, 1, grid, c, ζᶜᶜᶜ, mᵢᶜᶜᶜ, Δt)

    # Coefficient for substepping internal stress
    αᶜᶜᶜ = get_stepping_coefficients(i, j, 1, grid, substeps, c)
    αᶠᶠᶜ = ℑxyᴮᶠᶠᶜ(i, j, 1, grid, get_stepping_coefficients, substeps, c)

    @inbounds σ₁₁ₚ = (1 + γc) * σ₁₁[i, j, 1] - γc * ℑxyᴮᶜᶜᶜ(i, j, 1, grid, σ̂₁₁)
    @inbounds σ₂₂ₚ = (1 + γc) * σ₂₂[i, j, 1] - γc * ℑxyᴮᶜᶜᶜ(i, j, 1, grid, σ̂₂₂)
    @inbounds σ₁₂ₚ = (1 + γc) * σ₁₂[i, j, 1] - γc * ℑxyᴮᶠᶠᶜ(i, j, 1, grid, σ̂₁₂) 

    @inbounds σ̂₁₁ₚ = (1 + γc) * σ̂₁₁[i, j, 1] - γc * ℑxyᴮᶠᶠᶜ(i, j, 1, grid, σ₁₁)
    @inbounds σ̂₂₂ₚ = (1 + γc) * σ̂₂₂[i, j, 1] - γc * ℑxyᴮᶠᶠᶜ(i, j, 1, grid, σ₂₂)
    @inbounds σ̂₁₂ₚ = (1 + γc) * σ̂₁₂[i, j, 1] - γc * ℑxyᴮᶜᶜᶜ(i, j, 1, grid, σ₁₂) 

    @inbounds σ₁₁[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, (σ₁₁ᵖ⁺¹ - σ₁₁ₚ) / αᶜᶜᶜ, zero(grid))
    @inbounds σ₂₂[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, (σ₂₂ᵖ⁺¹ - σ₂₂ₚ) / αᶜᶜᶜ, zero(grid))
    @inbounds σ₁₂[i, j, 1] += ifelse(mᵢᶠᶠᶜ > 0, (σ₁₂ᵖ⁺¹ - σ₁₂ₚ) / αᶠᶠᶜ, zero(grid))

    @inbounds σ̂₁₁[i, j, 1] += ifelse(mᵢᶠᶠᶜ > 0, (σ̂₁₁ᵖ⁺¹ - σ̂₁₁ₚ) / αᶠᶠᶜ, zero(grid))
    @inbounds σ̂₂₂[i, j, 1] += ifelse(mᵢᶠᶠᶜ > 0, (σ̂₂₂ᵖ⁺¹ - σ̂₂₂ₚ) / αᶠᶠᶜ, zero(grid))
    @inbounds σ̂₁₂[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, (σ̂₁₂ᵖ⁺¹ - σ̂₁₂ₚ) / αᶜᶜᶜ, zero(grid))
end

#####
##### Internal stress divergence for the EVP model
#####

@inline function ∂ⱼ_σ₁ⱼᶠᶜᶜ(i, j, k, grid, ::ExplicitViscoPlasticRheology, fields) 
    ∂xσ₁₁ = ∂xᶠᶜᶜ_c(i, j, k, grid, fields.σ₁₁)
    ∂yσ₁₂ = ∂yᶠᶜᶜ_V(i, j, k, grid, fields.σ₁₂)

    return ∂xσ₁₁ + ∂yσ₁₂ 
end

@inline function ∂ⱼ_σ₁ⱼᶜᶠᶜ(i, j, k, grid, ::ExplicitViscoPlasticRheology, fields) 
    ∂xσ₁₁ = ∂xᶜᶠᶜ_U(i, j, k, grid, fields.σ̂₁₁)
    ∂yσ₁₂ = ∂yᶜᶠᶜ_c(i, j, k, grid, fields.σ̂₁₂)

    return ∂xσ₁₁ + ∂yσ₁₂ 
end

@inline function ∂ⱼ_σ₂ⱼᶜᶠᶜ(i, j, k, grid, ::ExplicitViscoPlasticRheology, fields) 
    ∂xσ₁₂ = ∂xᶜᶠᶜ_U(i, j, k, grid, fields.σ₁₂)
    ∂yσ₂₂ = ∂yᶜᶠᶜ_c(i, j, k, grid, fields.σ₂₂)

    return ∂xσ₁₂ + ∂yσ₂₂
end

@inline function ∂ⱼ_σ₂ⱼᶠᶜᶜ(i, j, k, grid, ::ExplicitViscoPlasticRheology, fields) 
    ∂xσ₁₂ = ∂xᶠᶜᶜ_c(i, j, k, grid, fields.σ̂₁₂)
    ∂yσ₂₂ = ∂yᶠᶜᶜ_V(i, j, k, grid, fields.σ̂₂₂)

    return ∂xσ₁₂ + ∂yσ₂₂
end

# To help convergence to the right velocities
@inline rheology_specific_forcing_xᶠᶜᶜ(i, j, k, grid, ::ExplicitViscoPlasticRheology, fields, uᵢ) = fields.uⁿ[i, j, k] - uᵢ[i, j, k]
@inline rheology_specific_forcing_yᶜᶠᶜ(i, j, k, grid, ::ExplicitViscoPlasticRheology, fields, vᵢ) = fields.vⁿ[i, j, k] - vᵢ[i, j, k]
@inline rheology_specific_forcing_xᶜᶠᶜ(i, j, k, grid, ::ExplicitViscoPlasticRheology, fields, ûᵢ) = fields.ûⁿ[i, j, k] - ûᵢ[i, j, k]
@inline rheology_specific_forcing_yᶠᶜᶜ(i, j, k, grid, ::ExplicitViscoPlasticRheology, fields, v̂ᵢ) = fields.v̂ⁿ[i, j, k] - v̂ᵢ[i, j, k]
