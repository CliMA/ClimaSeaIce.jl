struct BrittleBinghamMaxellRheology{FT, A}
    ice_ridging_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    ice_cohesion :: FT # cohesion
    undamaged_elastic_modulus :: FT # minimum plastic parameter (transitions to viscous behaviour)
    undamaged_viscous_relaxation_time :: FT # minimum number of substeps expressed as the dynamic coefficient
    ridging_ice_thickness :: FT # maximum number of substeps expressed as the dynamic coefficient
    poisson_ratio :: FT
    friction_coefficient :: FT
    maximum_compressive_strength :: FT
    healing_constant :: FT
    damage_parameter :: FT 
    damage_interpolation_scheme :: A # Interpolation of damage on face-face regions
end

function BrittleBinghamMaxellRheology(FT::DataType = Float64; 
                                      ice_ridging_strength = 1e4, 
                                      ice_compaction_hardening = 20, 
                                      ice_cohesion = 5.8e3,
                                      undamaged_elastic_modulus = 5.96e8,
                                      undamaged_viscous_relaxation_time = 1e7,
                                      ridging_ice_thickness = 1,
                                      poisson_ratio = 1 / 3,
                                      friction_coefficient = 0.7,
                                      maximum_compressive_strength = 2.9e7, 
                                      healing_constant = 26, # Ks
                                      damage_parameter = 5,
                                      damage_interpolation_scheme = Centered())

    return BrittleBinghamMaxellRheology(convert(FT, ice_ridging_strength), 
                                        convert(FT, ice_compaction_hardening), 
                                        convert(FT, ice_cohesion),
                                        convert(FT, undamaged_elastic_modulus),
                                        convert(FT, undamaged_viscous_relaxation_time),
                                        convert(FT, ridging_ice_thickness),
                                        convert(FT, poisson_ratio),
                                        convert(FT, friction_coefficient),
                                        convert(FT, maximum_compressive_strength),
                                        convert(FT, healing_constant),
                                        convert(FT, damage_parameter),
                                        damage_interpolation_scheme)
end

required_prognostic_tracers(::BrittleBinghamMaxellRheology, grid) = 
    (; d = Field{Center, Center, Nothing}(grid)) # damage tracer
    
function required_auxiliary_fields(::BrittleBinghamMaxellRheology, grid)
    
    # TODO: What about boundary conditions?
    P = Field{Center, Center, Nothing}(grid)
    E = Field{Center, Center, Nothing}(grid)
    λ = Field{Center, Center, Nothing}(grid)
    
    P̂ = Field{Face, Face, Nothing}(grid)
    Ê = Field{Face, Face, Nothing}(grid)
    λ̂ = Field{Face, Face, Nothing}(grid)
    
    σ₁₁ = Field{Center, Center, Nothing}(grid)
    σ₂₂ = Field{Center, Center, Nothing}(grid)
    σ₁₂ = Field{Face, Face, Nothing}(grid)

    σ̂₁₁ = Field{Face,   Face,   Nothing}(grid)
    σ̂₂₂ = Field{Face,   Face,   Nothing}(grid)
    σ̂₁₂ = Field{Center, Center, Nothing}(grid)

    d̂ = Field{Face, Face, Nothing}(grid)
    return (; σ₁₁, σ₂₂, σ₁₂, σ̂₁₁, σ̂₂₂, σ̂₁₂, d̂, P, E, λ, P̂, Ê, λ̂)
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
                                 Adapt.adapt(to, r.maximum_compressive_strength),
                                 Adapt.adapt(to, r.healing_constant),
                                 Adapt.adapt(to, r.damage_parameter),
                                 Adapt.adapt(to, r.damage_interpolation_scheme))

#####
##### Computation of the stresses
#####

function initialize_rheology!(model, rheology::BrittleBinghamMaxellRheology)
    h = model.ice_thickness
    ℵ = model.ice_concentration
    
    fields = model.dynamics.auxiliary_fields

    P★ = rheology.ice_ridging_strength
    E★ = rheology.undamaged_elastic_modulus
    λ★ = rheology.undamaged_viscous_relaxation_time
    h★ = rheology.ridging_ice_thickness
    α  = rheology.damage_parameter
    C  = rheology.ice_compaction_hardening
    
    # compute on the whole grid including halos
    parameters = KernelParameters(size(fields.P.data)[1:2], fields.P.data.offsets[1:2])
    launch!(architecture(model.grid), model.grid, parameters, _initialize_bbm_rheology!, fields, model.grid, P★, E★, λ★, h★, α, C, h, ℵ)
    
    return nothing
end

@kernel function _initialize_bbm_rheology!(fields, grid, P★, E★, λ★, h★, α, C, h, ℵ)
    i, j = @index(Global, NTuple)    
    @inbounds begin
        ecc = exp(- C * (1 - ℵ[i, j, 1])) 
        eff = exp(- C * (1 - ℑxyᶠᶠᵃ(i, j, 1, grid, ℵ)))
        
        # Center - Center fields
        fields.P[i, j, 1] = P★ * (h[i, j, 1] / h★)^(3/2) * ecc
        fields.E[i, j, 1] = E★ * ecc
        fields.λ[i, j, 1] = λ★ * ecc^(α - 1)

        # Face - Face fields
        fields.P̂[i, j, 1] = P★ * (ℑxyᶠᶠᵃ(i, j, 1, grid, h) / h★)^(3/2) * eff
        fields.Ê[i, j, 1] = E★ * eff
        fields.λ̂[i, j, 1] = λ★ * eff^(α - 1)
    end
end

# Specific compute stresses for the EVP rheology
function compute_stresses!(model, dynamics, rheology::BrittleBinghamMaxellRheology, Δt, Ns) 

    grid = model.grid
    arch = architecture(grid)

    ρᵢ   = model.ice_density
    d    = model.tracers.d
    u, v = model.velocities
    fields = dynamics.auxiliary_fields

    Nx, Ny, _ = size(grid)

    parameters = KernelParameters(-1:Nx+2, -1:Ny+2)

    # Pretty simple timestepping
    Δτ = Δt / Ns

    launch!(arch, grid, parameters, _advance_bbm_stresses!, fields, grid, rheology, d, u, v, Δτ)
    launch!(arch, grid, parameters, _mohr_colomb_correction!, fields, grid, rheology, d, ρᵢ, Δτ)

    return nothing
end

@inline σᴵᶜᶜᶜ(i, j, k, grid, fields) = @inbounds (fields.σ₁₁[i, j, k] + fields.σ₂₂[i, j, k]) / 2
@inline σᴵᶠᶠᶜ(i, j, k, grid, fields) = @inbounds (fields.σ̂₁₁[i, j, k] + fields.σ̂₂₂[i, j, k]) / 2

@inline function σᴵᴵᶜᶜᶜ(i, j, k, grid, fields) 
    σ₁₁ = @inbounds fields.σ₁₁[i, j, k]
    σ₂₂ = @inbounds fields.σ₂₂[i, j, k]
    σ₁₂ = @inbounds fields.σ̂₁₂[i, j, k]
    
    return sqrt((σ₁₁ - σ₂₂)^2 / 4 + σ₁₂^2)
end

@inline function σᴵᴵᶠᶠᶜ(i, j, k, grid, fields) 
    σ₁₁ = @inbounds fields.σ̂₁₁[i, j, k]
    σ₂₂ = @inbounds fields.σ̂₂₂[i, j, k]
    σ₁₂ = @inbounds fields.σ₁₂[i, j, k]
    
    return sqrt((σ₁₁ - σ₂₂)^2 / 4 + σ₁₂^2)
end

@kernel function _advance_bbm_stresses!(fields, grid, rheology, d, u, v, Δτ)
    i, j = @index(Global, NTuple)

    P = fields.P
    E = fields.E
    λ = fields.λ

    P̂ = fields.P̂
    Ê = fields.Ê
    λ̂ = fields.λ̂
    d̂ = fields.d̂

    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂

    σ̂₁₁ = fields.σ̂₁₁
    σ̂₂₂ = fields.σ̂₂₂
    σ̂₁₂ = fields.σ̂₁₂

    α = rheology.damage_parameter
    ν = rheology.poisson_ratio

    Pᶜᶜᶜ = @inbounds P[i, j, 1]
    Eᶜᶜᶜ = @inbounds E[i, j, 1] * (1 - d[i, j, 1])
    λᶜᶜᶜ = @inbounds λ[i, j, 1] * (1 - d[i, j, 1])^(α - 1)

    Pᶠᶠᶜ = @inbounds P̂[i, j, 1]
    Eᶠᶠᶜ = @inbounds Ê[i, j, 1] * (1 - d̂[i, j, 1])
    λᶠᶠᶜ = @inbounds λ̂[i, j, 1] * (1 - d̂[i, j, 1])^(α - 1)

    σIᶜᶜᶜ = σᴵᶜᶜᶜ(i, j, 1, grid, fields) 
    σIᶠᶠᶜ = σᴵᶠᶠᶜ(i, j, 1, grid, fields) 

    # Test which isotropic stress to use
    Pᶜᶜᶜ = ifelse(σIᶜᶜᶜ < - Pᶜᶜᶜ, Pᶜᶜᶜ / σIᶜᶜᶜ, 
           ifelse(σIᶜᶜᶜ > 0     , zero(grid), - one(grid)))

    # Test which isotropic stress to use
    Pᶠᶠᶜ = ifelse(σIᶠᶠᶜ < - Pᶠᶠᶜ, Pᶠᶠᶜ / σIᶠᶠᶜ, 
           ifelse(σIᶠᶠᶜ > 0     , zero(grid), - one(grid)))

    # Strain rates
    ϵ̇₁₁ = strain_rate_xx(i, j, 1, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, 1, grid, u, v) 
    ϵ̇₁₂ = strain_rate_xy(i, j, 1, grid, u, v)
    
    # Strain rates
    ϵ̇̂₁₁ = ℑxyᶠᶠᵃ(i, j, 1, grid, strain_rate_xx, u, v) 
    ϵ̇̂₂₂ = ℑxyᶠᶠᵃ(i, j, 1, grid, strain_rate_yy, u, v) 
    ϵ̇̂₁₂ = ℑxyᶜᶜᵃ(i, j, 1, grid, strain_rate_xy, u, v)
    
    Kϵ₁₁ = (ϵ̇₁₁ + ν * ϵ̇₂₂) / (1 - ν^2)
    Kϵ₂₂ = (ϵ̇₂₂ + ν * ϵ̇₁₁) / (1 - ν^2)
    Kϵ₁₂ =  (1 - ν) * ϵ̇₁₂  / (1 - ν^2)

    K̂ϵ₁₁ = (ϵ̇̂₁₁ + ν * ϵ̇̂₂₂) / (1 - ν^2)
    K̂ϵ₂₂ = (ϵ̇̂₂₂ + ν * ϵ̇̂₁₁) / (1 - ν^2)
    K̂ϵ₁₂ =  (1 - ν) * ϵ̇̂₁₂  / (1 - ν^2)

    # Implicit diagonal operator
    Ωᶜᶜᶜ = 1 / (1 + Δτ * (1 + Pᶜᶜᶜ) / λᶜᶜᶜ)
    Ωᶠᶠᶜ = 1 / (1 + Δτ * (1 + Pᶠᶠᶜ) / λᶠᶠᶜ)

    @inline σ₁₁[i, j, 1] = Ωᶜᶜᶜ * (σ₁₁[i, j, 1] + Δτ * Eᶜᶜᶜ * Kϵ₁₁)
    @inline σ₂₂[i, j, 1] = Ωᶜᶜᶜ * (σ₂₂[i, j, 1] + Δτ * Eᶜᶜᶜ * Kϵ₂₂)
    @inline σ₁₂[i, j, 1] = Ωᶠᶠᶜ * (σ₁₂[i, j, 1] + Δτ * Eᶠᶠᶜ * Kϵ₁₂)

    @inline σ̂₁₁[i, j, 1] = Ωᶠᶠᶜ * (σ̂₁₁[i, j, 1] + Δτ * Eᶠᶠᶜ * K̂ϵ₁₁)
    @inline σ̂₂₂[i, j, 1] = Ωᶠᶠᶜ * (σ̂₂₂[i, j, 1] + Δτ * Eᶠᶠᶜ * K̂ϵ₂₂)
    @inline σ̂₁₂[i, j, 1] = Ωᶜᶜᶜ * (σ̂₁₂[i, j, 1] + Δτ * Eᶜᶜᶜ * K̂ϵ₁₂)
end

@kernel function _mohr_colomb_correction!(fields, grid, rheology, d, ρᵢ, Δτ)
    i, j = @index(Global, NTuple)
    
    E = fields.E
    Ê = fields.Ê
    d̂ = fields.d̂
    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂

    σ̂₁₁ = fields.σ̂₁₁
    σ̂₂₂ = fields.σ̂₂₂
    σ̂₁₂ = fields.σ̂₁₂

    ν = rheology.poisson_ratio
    N = rheology.maximum_compressive_strength
    c = rheology.ice_cohesion
    μ = rheology.friction_coefficient

    Eᶜᶜᶜ = @inbounds E[i, j, 1] * (1 - d[i, j, 1])
    Eᶠᶠᶜ = @inbounds Ê[i, j, 1] * (1 - d̂[i, j, 1])

    # Principal stress invariant
    σIᶜᶜᶜ  =  σᴵᶜᶜᶜ(i, j, 1, grid, fields) 
    σIᶠᶠᶜ  =  σᴵᶠᶠᶜ(i, j, 1, grid, fields) 
    σIIᶜᶜᶜ = σᴵᴵᶜᶜᶜ(i, j, 1, grid, fields) 
    σIIᶠᶠᶜ = σᴵᴵᶠᶠᶜ(i, j, 1, grid, fields) 

    # critical damage computation
    dcᶜᶜᶜ = ifelse(σIᶜᶜᶜ > - N, c / (σIIᶜᶜᶜ + μ * σIᶜᶜᶜ), - N / σIᶜᶜᶜ)
    dcᶠᶠᶜ = ifelse(σIᶠᶠᶜ > - N, c / (σIIᶠᶠᶜ + μ * σIᶠᶠᶜ), - N / σIᶠᶠᶜ)

    # Relaxation time
    tdᶜᶜᶜ = @inbounds sqrt(2 * (1 + ν) * ρᵢ[i, j, 1] / Eᶜᶜᶜ * Azᶜᶜᶜ(i, j, 1, grid))
    tdᶠᶠᶜ = @inbounds sqrt(2 * (1 + ν) * ρᵢ[i, j, 1] / Eᶠᶠᶜ * Azᶠᶠᶜ(i, j, 1, grid))

    Gdᶜᶜ = @inbounds   (1 - dcᶜᶜᶜ) * (1 - d[i, j, 1]) * Δτ / tdᶜᶜᶜ
    Gdᶠᶠ = @inbounds   (1 - dcᶠᶠᶜ) * (1 - d̂[i, j, 1]) * Δτ / tdᶠᶠᶜ
    
    Gσ₁₁ = @inbounds - (1 - dcᶜᶜᶜ) * σ₁₁[i, j, 1]  * Δτ / tdᶜᶜᶜ
    Gσ₂₂ = @inbounds - (1 - dcᶜᶜᶜ) * σ₂₂[i, j, 1]  * Δτ / tdᶜᶜᶜ
    Gσ₁₂ = @inbounds - (1 - dcᶠᶠᶜ) * σ₁₂[i, j, 1]  * Δτ / tdᶠᶠᶜ
    
    Gσ̂₁₁ = @inbounds - (1 - dcᶠᶠᶜ) * σ̂₁₁[i, j, 1]  * Δτ / tdᶠᶠᶜ
    Gσ̂₂₂ = @inbounds - (1 - dcᶠᶠᶜ) * σ̂₂₂[i, j, 1]  * Δτ / tdᶠᶠᶜ
    Gσ̂₁₂ = @inbounds - (1 - dcᶜᶜᶜ) * σ̂₁₂[i, j, 1]  * Δτ / tdᶜᶜᶜ

    # Damage and stress updates
    @inbounds d[i, j, 1] += ifelse(0 ≤ dcᶜᶜᶜ ≤ 1, Gdᶜᶜ, zero(grid))
    @inbounds d̂[i, j, 1] += ifelse(0 ≤ dcᶠᶠᶜ ≤ 1, Gdᶠᶠ, zero(grid))

    @inbounds σ₁₁[i, j, 1] += ifelse(0 ≤ dcᶜᶜᶜ ≤ 1, Gσ₁₁, zero(grid))
    @inbounds σ₂₂[i, j, 1] += ifelse(0 ≤ dcᶜᶜᶜ ≤ 1, Gσ₂₂, zero(grid))
    @inbounds σ₁₂[i, j, 1] += ifelse(0 ≤ dcᶠᶠᶜ ≤ 1, Gσ₁₂, zero(grid))

    @inbounds σ̂₁₁[i, j, 1] += ifelse(0 ≤ dcᶠᶠᶜ ≤ 1, Gσ̂₁₁, zero(grid))
    @inbounds σ̂₂₂[i, j, 1] += ifelse(0 ≤ dcᶠᶠᶜ ≤ 1, Gσ̂₂₂, zero(grid))
    @inbounds σ̂₁₂[i, j, 1] += ifelse(0 ≤ dcᶜᶜᶜ ≤ 1, Gσ̂₁₂, zero(grid))

    # Clamp damage between 0 and 1
    @inbounds d[i, j, 1] = clamp(d[i, j, 1], zero(grid), 99999 * one(grid) / 100000)
    @inbounds d̂[i, j, 1] = clamp(d̂[i, j, 1], zero(grid), 99999 * one(grid) / 100000)
end

#####
##### Methods for the BBM rheology
#####

# In the BBM rheology, the stresses need to be vertically integrated
@inline hσ₁₁(i, j, k, grid, fields) = @inbounds fields.σ₁₁[i, j, k] * fields.h[i, j, k]
@inline hσ₂₂(i, j, k, grid, fields) = @inbounds fields.σ₂₂[i, j, k] * fields.h[i, j, k]
@inline hσ₁₂(i, j, k, grid, fields) = @inbounds fields.σ₁₂[i, j, k] * ℑxᶠᵃᵃ(i, j, k, grid, fields.h)

# Here we extend all the functions that a rheology model needs to support:
@inline function ∂ⱼ_σ₁ⱼ(i, j, k, grid, ::BrittleBinghamMaxellRheology, clock, fields) 
    ∂xσ₁₁ = δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶜᶜ, hσ₁₁, fields) / Azᶠᶜᶜ(i, j, k, grid)
    ∂yσ₁₂ = δyᵃᶜᵃ(i, j, k, grid, Δx_qᶠᶠᶜ, hσ₁₂, fields) / Azᶠᶜᶜ(i, j, k, grid)

    return ∂xσ₁₁ + ∂yσ₁₂
end

@inline function ∂ⱼ_σ₂ⱼ(i, j, k, grid, ::BrittleBinghamMaxellRheology, clock, fields) 
    ∂xσ₁₂ = δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶠᶜ, hσ₁₂, fields) / Azᶜᶠᶜ(i, j, k, grid)
    ∂yσ₂₂ = δyᵃᶠᵃ(i, j, k, grid, Δx_qᶜᶜᶜ, hσ₂₂, fields) / Azᶜᶠᶜ(i, j, k, grid)

    return ∂xσ₁₂ + ∂yσ₂₂
end