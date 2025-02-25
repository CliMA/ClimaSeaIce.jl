struct BrittleBinghamMaxwellRheology{FT}
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
end

function BrittleBinghamMaxwellRheology(FT::DataType = Float64; 
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
                                       damage_parameter = 5)

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
                                         convert(FT, damage_parameter))
end
    
function required_auxiliary_fields(::BrittleBinghamMaxwellRheology, grid)
    
    # TODO: What about boundary conditions?
    P = Field{Center, Center, Nothing}(grid)
    E = Field{Center, Center, Nothing}(grid)
    λ = Field{Center, Center, Nothing}(grid)
    d = Field{Center, Center, Nothing}(grid)

    σ₁₁ = Field{Center, Center, Nothing}(grid)
    σ₂₂ = Field{Center, Center, Nothing}(grid)
    σ₁₂ = Field{Center, Center, Nothing}(grid)

    return (; σ₁₁, σ₂₂, σ₁₂, P, E, λ, d)
end

# Extend the `adapt_structure` function for the `BrittleBinghamMaxwellRheology`
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
                                  Adapt.adapt(to, r.damage_parameter))

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
    u, v = model.velocities
    fields = dynamics.auxiliary_fields

    Nx, Ny, _ = size(grid)

    parameters = KernelParameters(-3:Nx+2, -3:Ny+2)

    # Pretty simple timestepping
    Δτ = Δt / Ns

    launch!(arch, grid, parameters, _advance_bbm_stresses!, fields, grid, rheology, u, v, Δτ)
    launch!(arch, grid, parameters, _mohr_colomb_correction!, fields, grid, rheology, ρᵢ, Δτ)

    return nothing
end

@inline function σᴵ(i, j, k, grid, fields) 
    σ₁₁ = @inbounds fields.σ₁₁[i, j, k]
    σ₂₂ = @inbounds fields.σ₂₂[i, j, k]

    return (σ₁₁ + σ₂₂) / 2
end

@inline function σᴵᴵ(i, j, k, grid, fields) 
    σ₁₁ = @inbounds fields.σ₁₁[i, j, k]
    σ₂₂ = @inbounds fields.σ₂₂[i, j, k]
    σ₁₂ = @inbounds fields.σ₁₂[i, j, k]
    
    return sqrt((σ₁₁ - σ₂₂)^2 / 4 + σ₁₂^2)
end

@inline function reconstruction_2d(i, j, k, grid, f, args...)
    fij = f(i,   j,   k, grid, args...)
    fmj = f(i-1, j,   k, grid, args...)
    fpj = f(i+1, j,   k, grid, args...)
    fim = f(i,   j-1, k, grid, args...)
    fip = f(i,   j+1, k, grid, args...)
    fmm = f(i-1, j-1, k, grid, args...)
    fmp = f(i-1, j+1, k, grid, args...)
    fpm = f(i+1, j-1, k, grid, args...)
    fpp = f(i+1, j+1, k, grid, args...)

    # remove NaNs
    isnanij = isnan(fij)
    isnanmj = isnan(fmj)
    isnanpj = isnan(fpj)
    isnanim = isnan(fim)
    isnanip = isnan(fip)
    isnanmm = isnan(fmm)
    isnanmp = isnan(fmp)
    isnanpm = isnan(fpm)
    isnanpp = isnan(fpp)

    fij = ifelse(isnanij, zero(grid), fij) / 4
    fmj = ifelse(isnanmj, zero(grid), fmj) / 8
    fpj = ifelse(isnanpj, zero(grid), fpj) / 8
    fim = ifelse(isnanim, zero(grid), fim) / 8
    fip = ifelse(isnanip, zero(grid), fip) / 8
    fmm = ifelse(isnanmm, zero(grid), fmm) / 16
    fmp = ifelse(isnanmp, zero(grid), fmp) / 16
    fpm = ifelse(isnanpm, zero(grid), fpm) / 16
    fpp = ifelse(isnanpp, zero(grid), fpp) / 16

    return (fij + fmj + fpj + fim + fip + fmm + fmp + fpm + fpp)
end

@kernel function _advance_bbm_stresses!(fields, grid, rheology, u, v, Δτ)
    i, j = @index(Global, NTuple)

    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂
    d   = fields.d

    α = rheology.damage_parameter
    ν = rheology.poisson_ratio
    
    P = @inbounds fields.P[i, j, 1]
    E = @inbounds fields.E[i, j, 1] * (1 - d[i, j, 1])
    λ = @inbounds fields.λ[i, j, 1] * (1 - d[i, j, 1])^(α - 1)

    # Test which isotropic stress to use
    P̃ = zero(grid) # ifelse(σI < - P, P / σI, 
                   # ifelse(σI > 0, zero(grid), - one(grid)))

    # Strain rates
    ϵ̇₁₁ = strain_rate_xx(i, j, 1, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, 1, grid, u, v) 
    ϵ̇₁₂ = ℑxyᶜᶜᵃ(i, j, 1, grid, strain_rate_xy, u, v)
    
    Kϵ₁₁ = (ϵ̇₁₁ + ν * ϵ̇₂₂) / (1 - ν^2)
    Kϵ₂₂ = (ϵ̇₂₂ + ν * ϵ̇₁₁) / (1 - ν^2)
    Kϵ₁₂ =  (1 - ν) * ϵ̇₁₂  / (1 - ν^2)

    # Implicit diagonal operator
    Ω = 1 / (1 + Δτ * (1 + P̃) / λ)

    @inline σ₁₁[i, j, 1] = Ω * (σ₁₁[i, j, 1] + Δτ * E * Kϵ₁₁)
    @inline σ₂₂[i, j, 1] = Ω * (σ₂₂[i, j, 1] + Δτ * E * Kϵ₂₂)
    @inline σ₁₂[i, j, 1] = Ω * (σ₁₂[i, j, 1] + Δτ * E * Kϵ₁₂)
end

@inline function critical_damage(i, j, k, grid, N, c, μ, fields)
     # Principal stress invariant
    σI  =  σᴵ(i, j, k, grid, fields) 
    σII = σᴵᴵ(i, j, k, grid, fields) 

    return ifelse(σI > - N, c / (σII + μ * σI), - N / σI)
end

@kernel function _mohr_colomb_correction!(fields, grid, rheology, ρᵢ, Δτ)
    i, j = @index(Global, NTuple)
    
    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂
    d   = fields.d

    ν = rheology.poisson_ratio
    N = rheology.maximum_compressive_strength
    c = rheology.ice_cohesion
    μ = rheology.friction_coefficient

    E = @inbounds fields.E[i, j, 1] * (1 - d[i, j, 1])
    
    # critical damage computation
    dcrit = reconstruction_2d(i, j, 1, grid, critical_damage, N, c, μ, fields)

    # Relaxation time
    td = @inbounds sqrt(2 * (1 + ν) * ρᵢ[i, j, 1] / E * Azᶜᶜᶜ(i, j, 1, grid))

    Gd   = @inbounds   (1 - dcrit) * (1 - d[i, j, 1]) * Δτ / td
    Gσ₁₁ = @inbounds - (1 - dcrit) *    σ₁₁[i, j, 1]  * Δτ / td
    Gσ₂₂ = @inbounds - (1 - dcrit) *    σ₂₂[i, j, 1]  * Δτ / td
    Gσ₁₂ = @inbounds - (1 - dcrit) *    σ₁₂[i, j, 1]  * Δτ / td

    # Damage and stress updates
    @inbounds   d[i, j, 1] += ifelse(0 ≤ dcrit ≤ 1, Gd  , zero(grid))
    @inbounds σ₁₁[i, j, 1] += ifelse(0 ≤ dcrit ≤ 1, Gσ₁₁, zero(grid))
    @inbounds σ₂₂[i, j, 1] += ifelse(0 ≤ dcrit ≤ 1, Gσ₂₂, zero(grid))
    @inbounds σ₁₂[i, j, 1] += ifelse(0 ≤ dcrit ≤ 1, Gσ₁₂, zero(grid))

    # Clamp damage between 0 and a value close to 1 (cannot do 1 because of the relaxation time)
    @inbounds d[i, j, 1] = clamp(d[i, j, 1], zero(grid), 99999 * one(grid) / 100000)
end

#####
##### Methods for the BBM rheology
#####

# In the BBM rheology, the stresses need to be vertically integrated
@inline hσ₁₁(i, j, k, grid, fields) = @inbounds fields.σ₁₁[i, j, k] * fields.h[i, j, k]
@inline hσ₂₂(i, j, k, grid, fields) = @inbounds fields.σ₂₂[i, j, k] * fields.h[i, j, k]
@inline hσ₁₂(i, j, k, grid, fields) = @inbounds fields.σ₁₂[i, j, k] * fields.h[i, j, k]

# using Oceananigans.Advection: _biased_interpolate_yᵃᶠᵃ, _biased_interpolate_xᶠᵃᵃ, LeftBias
@inline hσ₁₂ᶠᶠᶜ(i, j, k, grid, fields) = ℑxyᶠᶠᵃ(i, j, k, grid, hσ₁₂, fields)

@inline ice_stress_ux(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) = hσ₁₁(i, j, k, grid, fields)
@inline ice_stress_uy(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) = hσ₁₂ᶠᶠᶜ(i, j, k, grid, fields)
@inline ice_stress_vx(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) = hσ₁₂ᶠᶠᶜ(i, j, k, grid, fields)
@inline ice_stress_vy(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) = hσ₂₂(i, j, k, grid, fields)
