using Oceananigans.Grids: halo_size

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
                                         convert(FT, damage_parameter),
                                         nothing)
end

rheology_prognostic_tracers(::BrittleBinghamMaxwellRheology) =  (:d, :σ₁₁, :σ₁₂, :σ₂₂)
    
function rheology_auxiliary_fields(::BrittleBinghamMaxwellRheology, grid)
    
    # TODO: What about boundary conditions?
    P  = Field{Center, Center, Nothing}(grid)
    E  = Field{Center, Center, Nothing}(grid)
    λ  = Field{Center, Center, Nothing}(grid)
    σ¹ = Field{Center, Center, Nothing}(grid)
    σ² = Field{Center, Center, Nothing}(grid)
    return (; P, E, λ, σ¹, σ²)
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
    d  = model.tracers.d

    # compute on the whole grid including halos
    parameters = KernelParameters(size(fields.P.data)[1:2], fields.P.data.offsets[1:2])
    launch!(architecture(model.grid), model.grid, parameters, _initialize_bbm_rheology!, fields, d, model.grid, P★, E★, λ★, h★, α, C, h, ℵ)
    
    return nothing
end

@kernel function _initialize_bbm_rheology!(fields, d, grid, P★, E★, λ★, h★, α, C, h, ℵ)
    i, j = @index(Global, NTuple)    
    @inbounds begin
        ecc = exp(- C * (1 - ℵ[i, j, 1])) 
        # Center - Center fields
        fields.P[i, j, 1] = P★ * (h[i, j, 1] / h★)^(3/2) * ecc
        fields.E[i, j, 1] = E★ * ecc
        fields.λ[i, j, 1] = λ★ * ecc^(α - 1)
        # Heal the damage with a healing constant

        # Clamp the damage between 0 and a value close to 1
        dᵢ = d[i, j, 1]
        d[i, j, 1] = clamp(dᵢ, zero(dᵢ), 99999 * one(dᵢ) / 100000)
    end
end

# Specific compute stresses for the EVP rheology
function compute_stresses!(model, dynamics, rheology::BrittleBinghamMaxwellRheology, Δt, Ns) 

    grid = model.grid
    arch = architecture(grid)

    ρᵢ   = model.ice_density
    u, v = model.velocities
    fields = dynamics.auxiliary_fields
    tracers = model.tracers

    Nx, Ny, _ = size(grid)
    Hx, Hy, _ = halo_size(grid)

    parameters = KernelParameters(-Hx+2:Nx+Hx-1, -Hy+2:Ny+Hy-1)

    # Pretty simple timestepping
    Δτ = Δt / Ns

    launch!(arch, grid, parameters, _compute_stress_predictors!,   fields,  grid, rheology, tracers, u, v, ρᵢ, Δτ)
    launch!(arch, grid, parameters, _advance_stresses_and_damage!, tracers, grid, rheology, fields,  ρᵢ, Δτ)

    return nothing
end

@inline function critical_damage(i, j, k, grid, N, c, μ, fields)
    σ¹ = @inbounds fields.σ¹[i, j, k]
    σ² = @inbounds fields.σ²[i, j, k]
    dc = one(grid) - ifelse(σ¹ > - N, c / (σ² + μ * σ¹), - N / σ¹)
    dc = ifelse(isnan(dc), zero(grid), dc)
    return dc * (σ² > c - μ * σ¹)
end

@kernel function _compute_stress_predictors!(fields, grid, rheology, tracers, u, v, ρᵢ, Δτ)
    i, j = @index(Global, NTuple)

    α = rheology.damage_parameter
    ν = rheology.poisson_ratio
    
    ϵ̇₁₁ = strain_rate_xx(i, j, 1, grid, u, v) 
    ϵ̇₂₂ = strain_rate_yy(i, j, 1, grid, u, v) 
    ϵ̇₁₂ = strain_rate_xy(i, j, 1, grid, u, v)

    Kϵ₁₁ = (ϵ̇₁₁ + ν  * ϵ̇₂₂) / (1 - ν^2)
    Kϵ₂₂ = (ϵ̇₂₂ + ν  * ϵ̇₁₁) / (1 - ν^2)
    Kϵ₁₂ =   (1 - ν) * ϵ̇₁₂  / (1 - ν^2)

    dᵢ = @inbounds tracers.d[i, j, 1]
    P  = @inbounds fields.P[i, j, 1]
    E₀ = @inbounds fields.E[i, j, 1]
    λ₀ = @inbounds fields.λ[i, j, 1] 

    E  = E₀ * (1 - dᵢ)
    λ  = λ₀ * (1 - dᵢ)^(α - 1)

    # Test which isotropic stress to use
    P̃ = zero(grid) 

    # Implicit diagonal operator
    Ω = 1 / (1 + Δτ * (1 + P̃) / λ)

    @inbounds tracers.σ₁₁[i, j, 1] = Ω * (tracers.σ₁₁[i, j, 1] + Δτ * E * Kϵ₁₁)
    @inbounds tracers.σ₂₂[i, j, 1] = Ω * (tracers.σ₂₂[i, j, 1] + Δτ * E * Kϵ₂₂)
    @inbounds tracers.σ₁₂[i, j, 1] = Ω * (tracers.σ₁₂[i, j, 1] + Δτ * E * Kϵ₁₂)

    σ₁₁ = @inbounds tracers.σ₁₁[i, j, 1]
    σ₂₂ = @inbounds tracers.σ₂₂[i, j, 1]
    σ₁₂ = @inbounds tracers.σ₁₂[i, j, 1]

    # Principal stress components
    @inbounds fields.σ¹[i, j, 1] = (σ₁₁ + σ₂₂) / 2
    @inbounds fields.σ²[i, j, 1] = sqrt((σ₁₁ - σ₂₂)^2 / 4 + σ₁₂^2)
end

@kernel function _advance_stresses_and_damage!(tracers, grid, rheology, fields, ρᵢ, Δτ)
    i, j = @index(Global, NTuple)

    α = rheology.damage_parameter
    ν = rheology.poisson_ratio
    N = rheology.maximum_compressive_strength
    c = rheology.ice_cohesion
    μ = rheology.friction_coefficient

    ρ  = @inbounds ρᵢ[i, j, 1]
    dᵢ = @inbounds tracers.d[i, j, 1]
    E₀ = @inbounds fields.E[i, j, 1]

    # The relaxation time is time-dependent
    E  = E₀ * (1 - dᵢ)
    td = sqrt(2 * (1 + ν) * ρ / E * Azᶜᶜᶜ(i, j, 1, grid))   

    # damage tendency
    Gd = reconstruction_2d(i, j, 1, grid, critical_damage, N, c, μ, fields) / td
     
    # Damage update
    dᵢ += Gd * (1 - dᵢ) * Δτ 

    # Advance stresses and clamp damage
    @inbounds tracers.σ₁₁[i, j, 1] -= Gd * Δτ * tracers.σ₁₁[i, j, 1]
    @inbounds tracers.σ₂₂[i, j, 1] -= Gd * Δτ * tracers.σ₂₂[i, j, 1]
    @inbounds tracers.σ₁₂[i, j, 1] -= Gd * Δτ * tracers.σ₁₂[i, j, 1]
    @inbounds tracers.d[i, j, 1]    = clamp(dᵢ, zero(grid), 99999 * one(grid) / 100000)
end

#####
##### Methods for the BBM rheology
#####

# In the BBM rheology, the stresses need to be vertically integrated
@inline hσ₁₁(i, j, k, grid, fields) = @inbounds fields.σ₁₁[i, j, k] * fields.h[i, j, k]
@inline hσ₂₂(i, j, k, grid, fields) = @inbounds fields.σ₂₂[i, j, k] * fields.h[i, j, k]
@inline hσ₁₂(i, j, k, grid, fields) = @inbounds fields.σ₁₂[i, j, k] * fields.h[i, j, k]

# Here we extend all the functions that a rheology model needs to support:
@inline ice_stress_ux(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) = ℑyᵃᶠᵃ(i, j, 1, grid, hσ₁₁, fields)
@inline ice_stress_uy(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) = ℑxᶠᵃᵃ(i, j, 1, grid, hσ₁₂, fields)
@inline ice_stress_vx(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) = ℑyᵃᶠᵃ(i, j, 1, grid, hσ₁₂, fields)
@inline ice_stress_vy(i, j, k, grid, ::BrittleBinghamMaxwellRheology, clock, fields) = ℑxᶠᵃᵃ(i, j, 1, grid, hσ₂₂, fields)

# Another formulation of the critical Damage

# @inline function dc_perpendicular(i, j, k, grid, N, c, μ, fields)

#     σ¹ = @inbounds fields.σ¹[i, j, k]
#     σ² = @inbounds fields.σ²[i, j, k]

#     m = tand(90 - atand(μ))
#     q = σ² - m * σ¹

#     # # Move towards the yield curve in a perpendicular fashion
#     σᵪ¹ = (c - q) / (m + μ)
#     σᵪ² = m * σᵪ¹ + q

#     dc = one(grid) - sqrt(σᵪ¹^2 + σᵪ²^2) / sqrt(σ¹^2 + σ²^2)
#     dc = ifelse(isnan(dc), zero(grid), dc)

#     return dc * (σ² > c - μ * σ¹)
# end
#
# Dcrit needs to be reconstructed
#
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