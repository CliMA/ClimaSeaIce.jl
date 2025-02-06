struct BrittleBinghamMaxwellRheology{FT, A}
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
    interpolation_scheme :: A # Interpolation of cc variables to faces
end

function BrittleBinghamMaxwellRheology(FT::DataType = Float64; 
                                       ice_ridging_strength = 1e4, 
                                       ice_compaction_hardening = 20, 
                                       ice_cohesion = 5.7e3,
                                       undamaged_elastic_modulus = 5.96e8,
                                       undamaged_viscous_relaxation_time = 1e7,
                                       ridging_ice_thickness = 1,
                                       poisson_ratio = 1 / 3,
                                       friction_coefficient = 0.7,
                                       maximum_compressive_strength = 2.9e7, 
                                       healing_constant = 26, # Ks
                                       damage_parameter = 5,
                                       interpolation_scheme = Centered())

    return BrittleBinghamMaxwellRheology(convert(FT, ice_ridging_strength), 
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
                                         interpolation_scheme)
end

required_prognostic_tracers(::BrittleBinghamMaxwellRheology, grid) = 
    (; d = Field{Center, Center, Nothing}(grid)) # damage tracer
    
function required_auxiliary_fields(::BrittleBinghamMaxwellRheology, grid)
    
    # TODO: What about boundary conditions?
    P = Field{Center, Center, Nothing}(grid)
    E = Field{Center, Center, Nothing}(grid)
    λ = Field{Center, Center, Nothing}(grid)
    
    σ₁₁  = Field{Center, Center, Nothing}(grid)
    σ₂₂  = Field{Center, Center, Nothing}(grid)
    σ₁₂  = Field{Center, Center, Nothing}(grid)
    σ₁₁ₙ = Field{Center, Center, Nothing}(grid)
    σ₂₂ₙ = Field{Center, Center, Nothing}(grid)
    σ₁₂ₙ = Field{Center, Center, Nothing}(grid)
    uₙ   = Field{Face, Face, Center}(grid)
    vₙ   = Field{Face, Face, Center}(grid)
    dₙ   = Field{Center, Center, Nothing}(grid)
    σI   = Field{Center, Center, Nothing}(grid)
    σII  = Field{Center, Center, Nothing}(grid)
    σ₁₁ₒ = Field{Center, Center, Nothing}(grid)
    σ₂₂ₒ = Field{Center, Center, Nothing}(grid)
    σ₁₂ₒ = Field{Center, Center, Nothing}(grid)
    dₒ   = Field{Center, Center, Nothing}(grid)
    dcr  = Field{Center, Center, Nothing}(grid)
    td   = Field{Center, Center, Nothing}(grid)
    dadd = Field{Center, Center, Nothing}(grid)
    return (; σ₁₁, σ₂₂, σ₁₂, σ₁₁ₙ, σ₂₂ₙ, σ₁₂ₙ, dₙ, P, E, λ, uₙ, vₙ, σ₁₁ₒ, σ₂₂ₒ, σ₁₂ₒ, dₒ, dcr, σI, σII, td, dadd)
end

# Extend the `adapt_structure` function for the ElastoViscoPlasticRheology
Adapt.adapt_structure(to, r::BrittleBinghamMaxwellRheology) = 
    BrittleBinghamMaxwellRheology(Adapt.adapt(to, r.ice_ridging_strength), 
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
                                 Adapt.adapt(to, r.interpolation_scheme))

#####
##### Computation of the stresses
#####

function initialize_rheology!(model, rheology::BrittleBinghamMaxwellRheology)
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
        # Center - Center fields
        fields.P[i, j, 1] = P★ * (h[i, j, 1] / h★)^(3/2) * ecc
        fields.E[i, j, 1] = E★ * ecc
        fields.λ[i, j, 1] = λ★ * ecc^(α - 1)
    end
end

# Specific compute stresses for the EVP rheology
function compute_stresses!(model, dynamics, rheology::BrittleBinghamMaxwellRheology, Δt, Ns) 

    grid = model.grid
    arch = architecture(grid)

    ρᵢ   = model.ice_density
    d    = model.tracers.d
    u, v = model.velocities
    fields = dynamics.auxiliary_fields

    Nx, Ny, _ = size(grid)

    parameters = KernelParameters(-6:Nx+7, -6:Ny+7)

    # Pretty simple timestepping
    Δτ = Δt / Ns

    launch!(arch, grid, parameters, _advance_bbm_stresses!, fields, grid, rheology, d, u, v, ρᵢ, Δτ)

    return nothing
end

@inline strain_rate_xx(i, j, k, grid, u, v) = ℑyᵃᶜᵃ(i, j, k, grid, δxᶜᵃᵃ, Δy_qᶠᶠᶜ, u) / Azᶜᶜᶜ(i, j, k, grid)
@inline strain_rate_yy(i, j, k, grid, u, v) = ℑxᶜᵃᵃ(i, j, k, grid, δyᵃᶜᵃ, Δx_qᶠᶠᶜ, v) / Azᶜᶜᶜ(i, j, k, grid)

@inline strain_rate_xy(i, j, k, grid, u, v) = 
        (ℑyᵃᶜᵃ(i, j, k, grid, δxᶜᵃᵃ, Δy_qᶠᶠᶜ, v) + 
         ℑxᶜᵃᵃ(i, j, k, grid, δyᵃᶜᵃ, Δx_qᶠᶠᶜ, u)) / Azᶜᶜᶜ(i, j, k, grid) / 2

@inline  σᴵ(i, j, k, grid, fields) = @inbounds (fields.σ₁₁[i, j, k] + fields.σ₂₂[i, j, k]) / 2
@inline σᴵᴵ(i, j, k, grid, fields) = @inbounds sqrt((fields.σ₁₁[i, j, k] - fields.σ₂₂[i, j, k])^2 / 4 + fields.σ₁₂[i, j, k]^2)

@inline function dcritical(i, j, k, grid, N, c, μ, fields)
    σI  = σᴵ(i, j, k, grid, fields)
    σII = σᴵᴵ(i, j, k, grid, fields)
    return ifelse(σI > - N, c / (σII + μ * σI), - N / σI)
end

@kernel function _advance_bbm_stresses2!(fields, grid, rheology, d, u, v, ρᵢ, Δτ)
    i, j = @index(Global, NTuple)

    α = rheology.damage_parameter
    ν = rheology.poisson_ratio
    N = rheology.maximum_compressive_strength
    c = rheology.ice_cohesion
    μ = rheology.friction_coefficient
    
    ϵ̇₁₁ = strain_rate_xx(i, j, 1, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, 1, grid, u, v) 
    ϵ̇₁₂ = strain_rate_xy(i, j, 1, grid, u, v)

    Kϵ₁₁ = (ϵ̇₁₁ + ν  * ϵ̇₂₂) / (1 - ν^2)
    Kϵ₂₂ = (ϵ̇₂₂ + ν  * ϵ̇₁₁) / (1 - ν^2)
    Kϵ₁₂ =   (1 - ν) * ϵ̇₁₂  / (1 - ν^2)

    dᵢ = @inbounds d[i, j, 1]
    ρ  = @inbounds ρᵢ[i, j, 1]
    P  = @inbounds fields.P[i, j, 1]
    σ₁ = @inbounds fields.σ₁₁[i, j, 1]
    σ₂ = @inbounds fields.σ₂₂[i, j, 1]
    σ₃ = @inbounds fields.σ₁₂[i, j, 1]
    E₀ = @inbounds fields.E[i, j, 1]
    λ₀ = @inbounds fields.λ[i, j, 1] 

    σI  = (σ₁ + σ₂) / 2
    σII = sqrt((σ₁ - σ₂)^2 / 4 + σ₃^2)

    σ₀₁ = @inbounds fields.σ₁₁ₙ[i, j, 1]
    σ₀₂ = @inbounds fields.σ₂₂ₙ[i, j, 1]
    σ₀₃ = @inbounds fields.σ₁₂ₙ[i, j, 1]
    d₀  = @inbounds fields.dₙ[i, j, 1]

    E  = E₀ * (1 - dᵢ)
    λ  = λ₀ * (1 - dᵢ)^(α - 1)

    # Test which isotropic stress to use
    P̃ = clamp(P / σI, -one(grid), zero(grid))

    # Implicit diagonal operator
    Ω = 1 / (1 + Δτ * (1 + P̃) / λ)

    σ₁ =  Ω * (σ₀₁ + Δτ * E * Kϵ₁₁)
    σ₂ =  Ω * (σ₀₂ + Δτ * E * Kϵ₂₂)
    σ₃ =  Ω * (σ₀₃ + Δτ * E * Kϵ₁₂)

    σI  = (σ₁ + σ₂) / 2
    σII = sqrt((σ₁ - σ₂)^2 / 4 + σ₃^2)

    dcrit = ifelse(σI > - N, c / (σII + μ * σI), - N / σI) 

    # Relaxation time
    td = sqrt(2 * (1 + ν) * ρ / E * Azᶜᶜᶜ(i, j, 1, grid))        
    Gd = (1 - dcrit) * Δτ / td

    # # Damage and stress updates
    dᵢ = d₀ + ifelse(0 ≤ dcrit ≤ 1, Gd * (1 - dᵢ), zero(grid))
    σ₁ -= ifelse(0 ≤ dcrit ≤ 1, Gd * σ₁, zero(grid))
    σ₂ -= ifelse(0 ≤ dcrit ≤ 1, Gd * σ₂, zero(grid))
    σ₃ -= ifelse(0 ≤ dcrit ≤ 1, Gd * σ₃, zero(grid))

    # Clamp damage between 0 and a value close to 1 (cannot do 1 because of the relaxation time)
    dᵢ = clamp(dᵢ, zero(grid), 99999 * one(grid) / 100000)

    @inbounds  fields.dcr[i, j, 1] = dcrit
    @inbounds   fields.σI[i, j, 1] = σI
    @inbounds  fields.σII[i, j, 1] = σII
    @inbounds fields.dadd[i, j, 1] = 0 ≤ dcrit ≤ 1
    @inbounds   fields.td[i, j, 1] = td
    @inbounds  fields.σ₁₁[i, j, 1] = σ₁
    @inbounds  fields.σ₂₂[i, j, 1] = σ₂
    @inbounds  fields.σ₁₂[i, j, 1] = σ₃
    @inbounds           d[i, j, 1] = dᵢ
end

@kernel function _advance_bbm_stresses!(fields, grid, rheology, d, u, v, ρᵢ, Δτ)
    i, j = @index(Global, NTuple)

    α = rheology.damage_parameter
    ν = rheology.poisson_ratio
    I = rheology.interpolation_scheme

    ν = rheology.poisson_ratio
    N = rheology.maximum_compressive_strength
    c = rheology.ice_cohesion
    μ = rheology.friction_coefficient
    
    ϵ̇₁₁ = strain_rate_xx(i, j, 1, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, 1, grid, u, v) 
    ϵ̇₁₂ = strain_rate_xy(i, j, 1, grid, u, v)

    Kϵ₁₁ = (ϵ̇₁₁ + ν  * ϵ̇₂₂) / (1 - ν^2)
    Kϵ₂₂ = (ϵ̇₂₂ + ν  * ϵ̇₁₁) / (1 - ν^2)
    Kϵ₁₂ =   (1 - ν) * ϵ̇₁₂  / (1 - ν^2)

    dᵢ = @inbounds d[i, j, 1]
    ρ  = @inbounds ρᵢ[i, j, 1]
    P  = @inbounds fields.P[i, j, 1]
    σ₁ = @inbounds fields.σ₁₁[i, j, 1]
    σ₂ = @inbounds fields.σ₂₂[i, j, 1]
    σ₃ = @inbounds fields.σ₁₂[i, j, 1]
    E₀ = @inbounds fields.E[i, j, 1]
    λ₀ = @inbounds fields.λ[i, j, 1] 
    
    dcrit = 1.0

    σI  = (σ₁ + σ₂) / 2
    σII = sqrt((σ₁ - σ₂)^2 / 4 + σ₃^2)
    td  = Δτ

    σ₀₁ = σ₁
    σ₀₂ = σ₂
    σ₀₃ = σ₃

    E  = E₀ * (1 - dᵢ)
    λ  = λ₀ * (1 - dᵢ)^(α - 1)

    # Test which isotropic stress to use
    P̃ = clamp(P / σI, -one(grid), zero(grid))

    # Implicit diagonal operator
    Ω = 1 / (1 + Δτ * (1 + P̃) / λ)

    σ₁ =  Ω * (σ₀₁ + Δτ * E * Kϵ₁₁)
    σ₂ =  Ω * (σ₀₂ + Δτ * E * Kϵ₂₂)
    σ₃ =  Ω * (σ₀₃ + Δτ * E * Kϵ₁₂)

    σI  = (σ₁ + σ₂) / 2
    σII = sqrt((σ₁ - σ₂)^2 / 4 + σ₃^2)

    dcrit = ifelse(σI > - N, c / (σII + μ * σI), - N / σI) 

    # Relaxation time
    td = sqrt(2 * (1 + ν) * ρ / E * Azᶜᶜᶜ(i, j, 1, grid))        
    Gd = (1 - dcrit) * Δτ / td

    # # Damage and stress updates
    dᵢ = d[i, j, 1] + ifelse(0 ≤ dcrit ≤ 1, Gd * (1 - dᵢ), zero(grid))
    σ₁ -= ifelse(0 ≤ dcrit ≤ 1, Gd * σ₁, zero(grid))
    σ₂ -= ifelse(0 ≤ dcrit ≤ 1, Gd * σ₂, zero(grid))
    σ₃ -= ifelse(0 ≤ dcrit ≤ 1, Gd * σ₃, zero(grid))

    # Clamp damage between 0 and a value close to 1 (cannot do 1 because of the relaxation time)
        dᵢ = clamp(dᵢ, zero(grid), 99999 * one(grid) / 100000)

    @inbounds  fields.dcr[i, j, 1] = dcrit
    @inbounds   fields.σI[i, j, 1] = σI
    @inbounds  fields.σII[i, j, 1] = σII
    @inbounds fields.dadd[i, j, 1] = 0 ≤ dcrit ≤ 1
    @inbounds   fields.td[i, j, 1] = td
    @inbounds  fields.σ₁₁[i, j, 1] = σ₁
    @inbounds  fields.σ₂₂[i, j, 1] = σ₂
    @inbounds  fields.σ₁₂[i, j, 1] = σ₃
    @inbounds           d[i, j, 1] = dᵢ
end

#####
##### Methods for the BBM rheology
#####

# In the BBM rheology, the stresses need to be vertically integrated
@inline hσ₁₁(i, j, k, grid, fields) = @inbounds fields.σ₁₁[i, j, k] * fields.h[i, j, k]
@inline hσ₂₂(i, j, k, grid, fields) = @inbounds fields.σ₂₂[i, j, k] * fields.h[i, j, k]
@inline hσ₁₂(i, j, k, grid, fields) = @inbounds fields.σ₁₂[i, j, k] * fields.h[i, j, k]

# Here we extend all the functions that a rheology model needs to support:
@inline function ∂ⱼ_σ₁ⱼ(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) 

    ∂xσ₁₁  = ℑyᵃᶠᵃ(i, j, k, grid, δxᶠᵃᵃ, Δy_qᶜᶜᶜ, hσ₁₁, fields) / Azᶠᶠᶜ(i, j, k, grid)
    ∂yσ₁₂  = ℑxᶠᵃᵃ(i, j, k, grid, δyᵃᶠᵃ, Δx_qᶜᶜᶜ, hσ₁₂, fields) / Azᶠᶠᶜ(i, j, k, grid)

    return ∂xσ₁₁ + ∂yσ₁₂
end

@inline function ∂ⱼ_σ₂ⱼ(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) 
    
    ∂xσ₁₂  = ℑyᵃᶠᵃ(i, j, k, grid, δxᶠᵃᵃ, Δy_qᶜᶜᶜ, hσ₁₂, fields) / Azᶠᶠᶜ(i, j, k, grid)
    ∂yσ₂₂  = ℑxᶠᵃᵃ(i, j, k, grid, δyᵃᶠᵃ, Δx_qᶜᶜᶜ, hσ₂₂, fields) / Azᶠᶠᶜ(i, j, k, grid)

    return ∂xσ₁₂ + ∂yσ₂₂
end