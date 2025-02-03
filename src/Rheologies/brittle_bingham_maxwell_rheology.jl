struct BrittleBinghamMaxellRheology{FT, A}
    ice_ridging_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    ice_cohesion :: FT # cohesion
    undamaged_elastic_modulus :: FT # minimum plastic parameter (transitions to viscous behaviour)
    undamaged_viscous_relaxation_time :: FT # minimum number of substeps expressed as the dynamic coefficient
    poisson_ratio :: FT
    friction_coefficient :: FT
    damage_parameter :: FT 
    ridging_ice_thickness :: FT # maximum number of substeps expressed as the dynamic coefficient
    damage_propagation_timescale :: FT
    damage_interpolation_scheme :: A # Interpolation of damage on face-face regions
end

function BrittleBinghamMaxellRheology(FT::DataType = Float64; 
                                    ice_ridging_strength = 1e4, 
                                    ice_compaction_hardening = 20, 
                                    ice_cohesion = 5.8e3,
                                    undamaged_elastic_modulus = 5.96e8,
                                    undamaged_viscous_relaxation_time = 1e7,
                                    ridging_ice_thickness = 1,
                                    poisson_ratio = 1/3,
                                    friction_coefficient = 0.5,
                                    damage_parameter = 5,
                                    damage_propagation_timescale = 1e-3,
                                    damage_interpolation_scheme = Centered())

    return BrittleBinghamMaxellRheology(convert(FT, ice_ridging_strength), 
                                        convert(FT, ice_compaction_hardening), 
                                        convert(FT, ice_cohesion),
                                        convert(FT, undamaged_elastic_modulus),
                                        convert(FT, undamaged_viscous_relaxation_time),
                                        convert(FT, ridging_ice_thickness),
                                        convert(FT, poisson_ratio),
                                        convert(FT, friction_coefficient),
                                        convert(FT, damage_parameter),
                                        convert(FT, damage_propagation_timescale),
                                        damage_interpolation_scheme)
end

required_prognostic_tracers(::BrittleBinghamMaxellRheology, grid) = 
    (; d = Field{Center, Center, Nothing}(grid)) # damage tracer
    
function required_auxiliary_fields(::BrittleBinghamMaxellRheology, grid)
    
    # TODO: What about boundary conditions?
    P  = Field{Center, Center, Nothing}(grid)
    E  = Field{Center, Center, Nothing}(grid)
    λ  = Field{Center, Center, Nothing}(grid)

    σ₁₁ = Field{Center, Center, Nothing}(grid)
    σ₂₂ = Field{Center, Center, Nothing}(grid)
    σ₁₂ = Field{Face, Face, Nothing}(grid)

    return (; σ₁₁, σ₂₂, σ₁₂, P, E, λ)
end

# Extend the `adapt_structure` function for the ElastoViscoPlasticRheology
Adapt.adapt_structure(to, r::BrittleBinghamMaxellRheology) = 
    BrittleBinghamMaxellRheology(Adapt.adapt(to, r.ice_ridging_strength), 
                                 Adapt.adapt(to, r.ice_compaction_hardening), 
                                 Adapt.adapt(to, r.ice_cohesion),
                                 Adapt.adapt(to, r.undamaged_elastic_modulus),
                                 Adapt.adapt(to, r.undamaged_viscous_relaxation_time),
                                 Adapt.adapt(to, r.ridging_ice_thickness),
                                 Adapt.adapt(to, r.poisson_ratio),
                                 Adapt.adapt(to, r.friction_coefficient),
                                 Adapt.adapt(to, r.damage_parameter),
                                 Adapt.adapt(to, r.damage_propagation_timescale),
                                 Adapt.adapt(to, r.damage_interpolation_scheme))

#####
##### Computation of the stresses
#####

"""
    initialize_rheology!(model, rheology::BrittleBinghamMaxellRheology)

Initialize the brittle-bingham-maxwell rheology.
In this step we calculate the ice strength given the ice mass (thickness and concentration).
"""
function initialize_rheology!(model, rheology::BrittleBinghamMaxellRheology)
    h = model.ice_thickness
    ℵ = model.ice_concentration

    P★ = rheology.ice_ridging_strength
    E★ = rheology.undamaged_elastic_modulus
    λ★ = rheology.undamaged_viscous_relaxation_time
    h★ = rheology.ridging_ice_thickness
    α  = rheology.damage_parameter
    C  = rheology.ice_compaction_hardening
    
    fields = model.dynamics.auxiliary_fields

    # compute on the whole grid including halos
    parameters = KernelParameters(size(fields.P.data)[1:2], fields.P.data.offsets[1:2])
    launch!(architecture(model.grid), model.grid, parameters, _initialize_bbm_rhology!, fields, P★, E★, λ★, h★, α, C, h, ℵ)
    
    return nothing
end

@kernel function _initialize_bbm_rhology!(fields, P★, E★, λ★, h★, α, C, h, ℵ)
    i, j = @index(Global, NTuple)    
    @inbounds exponent = exp(- C * (1 - ℵ[i, j, k])) 
    @inbounds fields.P[i, j, 1] = P★ * (h[i, j, k] / h★)^(3/2) * exponent
    @inbounds fields.E[i, j, 1] = E★ * exponent
    @inbounds fields.λ[i, j, 1] = λ★ * exponent^(α - 1)
end

# Specific compute stresses for the EVP rheology
function compute_stresses!(model, dynamics, rheology::ElastoViscoPlasticRheology, Δt, Ns) 

    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ρᵢ = model.ice_density
    ℵ  = model.ice_concentration
    d  = model.tracers.d
    u, v = model.velocities
    fields = dynamics.auxiliary_fields

    Nx, Ny, _ = size(grid)

    parameters = KernelParameters(-1:Nx+2, -1:Ny+2)

    # Pretty simple timestepping
    Δτ = Δt / Ns

    launch!(arch, grid, parameters, _advance_bbm_stresses!, fields, grid, rheology, d, u, v, Δτ)
    launch!(arch, grid, parameters, _advance_damage_correct_stresses!, fields, grid, rheology, d, u, v, h, ℵ, ρᵢ, Δt)

    return nothing
end
 
# Very simple reconstruction (4-point average)
@inline reconstruct_on_nodes(i, j, k, grid, scheme, d) = ℑxyᶠᶠᵃ(i, j, k, grid, d)

@kernel function _advance_bbm_stresses!(fields, grid, rheology, d, u, v, Δτ)
    i, j = @index(Global, NTuple)

    P = fields.P
    E = fields.E
    λ = fields.λ

    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂

    α = rheology.damage_parameter
    ν = rheology.poisson_ratio

    Pᶜᶜᶜ = @inbounds P[i, j, 1]
    Eᶜᶜᶜ = @inbounds E[i, j, 1] * (1 - d[i, j, 1])
    λᶜᶜᶜ = @inbounds λ[i, j, 1] * (1 - d[i, j, 1])^(α - 1)

    dᶠᶠᶜ = reconstruct_on_nodes(i, j, 1, grid, rheology.damage_interpolation_scheme, d)

    Pᶠᶠᶜ = @inbounds ℑxyᶠᶠᵃ(i, j, 1, grid, P)
    Eᶠᶠᶜ = @inbounds ℑxyᶠᶠᵃ(i, j, 1, grid, E) * (1 - dᶠᶠᶜ)
    λᶠᶠᶜ = @inbounds ℑxyᶠᶠᵃ(i, j, 1, grid, λ) * (1 - dᶠᶠᶜ)^(α - 1)

    # Strain rates
    ϵ̇₁₁ = strain_rate_xx(i, j, 1, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, 1, grid, u, v) 
    ϵ̇₁₂ = strain_rate_xy(i, j, 1, grid, u, v)
    
    Kϵ₁₁ = (ϵ̇₁₁ + ν * ϵ̇₂₂) / (1 - ν^2)
    Kϵ₂₂ = (ϵ̇₂₂ + ν * ϵ̇₁₁) / (1 - ν^2)
    Kϵ₁₂ = (1 - ν) * ϵ̇₁₂

    @inline σ₁₁[i, j, 1] += Δτ * (Eᶜᶜᶜ * Kϵ₁₁ - σ₁₁[i, j, 1] / λᶜᶜᶜ * (1 + Pᶜᶜᶜ))
    @inline σ₂₂[i, j, 1] += Δτ * (Eᶜᶜᶜ * Kϵ₂₂ - σ₂₂[i, j, 1] / λᶜᶜᶜ * (1 + Pᶜᶜᶜ))
    @inline σ₁₂[i, j, 1] += Δτ * (Eᶠᶠᶜ * Kϵ₁₂ - σ₁₂[i, j, 1] / λᶠᶠᶜ * (1 + Pᶠᶠᶜ))
end

# Compute the visco-plastic stresses for a slab sea ice model.
# The function updates the internal stress variables `σ₁₁`, `σ₂₂`, and `σ₁₂` in the `rheology` object
# following the αEVP formulation of Kimmritz et al (2016).
@kernel function _compute_evp_stresses!(fields, grid, rheology, u, v, h, ℵ, ρᵢ, Δt)
    i, j = @index(Global, NTuple)

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.minimum_plastic_stress
    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂
    α   = fields.α

    # Strain rates
    ϵ̇₁₁ = strain_rate_xx(i, j, 1, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, 1, grid, u, v) 
    ϵ̇₁₂ = strain_rate_xy(i, j, 1, grid, u, v)

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
    γ²ᶜᶜᶜ = ifelse(isnan(γ²ᶜᶜᶜ), rheology.max_substeps^2, γ²ᶜᶜᶜ) # In case both ζᶜᶜᶜ and mᵢᶜᶜᶜ are zero
    γᶜᶜᶜ  = clamp(sqrt(γ²ᶜᶜᶜ), rheology.min_substeps, rheology.max_substeps)

    γ²ᶠᶠᶜ = ζᶠᶠᶜ * π^2 * Δt / mᵢᶠᶠᶜ / Azᶠᶠᶜ(i, j, 1, grid)
    γ²ᶠᶠᶜ = ifelse(isnan(γ²ᶠᶠᶜ), rheology.max_substeps^2, γ²ᶠᶠᶜ) # In case both ζᶠᶠᶜ and mᵢᶠᶠᶜ are zero
    γᶠᶠᶜ  = clamp(sqrt(γ²ᶠᶠᶜ), rheology.min_substeps, rheology.max_substeps)

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
##### Methods for the BBM rheology
#####

# Here we extend all the functions that a rheology model needs to support:
@inline function ∂ⱼ_σ₁ⱼ(i, j, k, grid, ::BrittleBinghamMaxellRheology, clock, fields) 
    ∂xσ₁₁ = δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶜᶜ, fields.σ₁₁) / Azᶠᶜᶜ(i, j, k, grid)
    ∂yσ₁₂ = δyᵃᶜᵃ(i, j, k, grid, Δx_qᶠᶠᶜ, fields.σ₁₂) / Azᶠᶜᶜ(i, j, k, grid)

    return ∂xσ₁₁ + ∂yσ₁₂
end

@inline function ∂ⱼ_σ₂ⱼ(i, j, k, grid, ::BrittleBinghamMaxellRheology, clock, fields) 
    ∂xσ₁₂ = δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶠᶜ, fields.σ₁₁) / Azᶜᶠᶜ(i, j, k, grid)
    ∂yσ₂₂ = δyᵃᶠᵃ(i, j, k, grid, Δx_qᶜᶜᶜ, fields.σ₁₂) / Azᶜᶠᶜ(i, j, k, grid)

    return ∂xσ₁₂ + ∂yσ₂₂
end