
struct BrittleBinghamMaxwellRheology{S1, S2, S3, E, H, D, P, FT} <: AbstractRheology
    σ₁₁ :: S1 # internal stress xx
    σ₂₂ :: S2 # internal stress yy
    σ₁₂ :: S3 # internal stress xy
    young_modulus :: E  # Young modulus
    bulk_viscosity :: H  # Bulk viscosity
    damage :: D # field containing ice damage
    ice_strength :: P # field containing the precomputed ice strength
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    undamaged_young_modulus  :: FT
    undamaged_bulk_viscosity :: FT
    bulk_viscosity_exponent :: FT
end

function BrittleBinghamMaxwellRheology(grid::AbstractGrid; 
                                       ice_compressive_strength = 27500, 
                                       ice_compaction_hardening = 20,
                                       undamaged_young_modulus = 100,
                                       undamaged_bulk_viscosity = 1e5,
                                       bulk_viscosity_exponent = 0.1)
    σ₁₁ = CenterField(grid)
    σ₂₂ = CenterField(grid)
    σ₁₂ = Field{Face, Face, Center}(grid)
    
    D  = CenterField(grid)
    E  = CenterField(grid)
    H  = CenterField(grid)
    P  = CenterField(grid)

    FT = eltype(grid)
    return BrittleBinghamMaxwellRheology(σ₁₁, σ₂₂, σ₁₂, E, H, D, P, 
                                         convert(FT, ice_compressive_strength), 
                                         convert(FT, ice_compaction_hardening),
                                         convert(FT, undamaged_young_modulus), 
                                         convert(FT, undamaged_bulk_viscosity), 
                                         convert(FT, bulk_viscosity_exponent))
end

# Extend the `adapt_structure` function for the BrittleBinghamMaxwellRheology
Adapt.adapt_structure(to, r::BrittleBinghamMaxwellRheology) = 
    BrittleBinghamMaxwellRheology(Adapt.adapt(to, r.σ₁₁),
                                  Adapt.adapt(to, r.σ₂₂),
                                  Adapt.adapt(to, r.σ₁₂),
                                  Adapt.adapt(to, r.young_modulus),
                                  Adapt.adapt(to, r.bulk_viscosity),
                                  Adapt.adapt(to, r.damage),
                                  Adapt.adapt(to, r.ice_strength),
                                  Adapt.adapt(to, r.ice_compressive_strength),
                                  Adapt.adapt(to, r.ice_compaction_hardening),
                                  Adapt.adapt(to, r.undamaged_young_modulus),
                                  Adapt.adapt(to, r.undamaged_bulk_viscosity),
                                  Adapt.adapt(to, r.bulk_viscosity_exponent))



"""
    initialize_rheology!(model, rheology::BrittleBinghamMaxwellRheology)

Initialize the elasto-visco-plastic rheology.
In this step we calculate the ice strength given the ice mass (thickness and concentration).
"""
function initialize_rheology!(model, rheology::BrittleBinghamMaxwellRheology)
    h = model.ice_thickness
    ℵ = model.ice_concentration

    P  = rheology.ice_strength
    E  = rheology.young_modulus
    η  = rheology.bulk_viscosity

    P★ = rheology.ice_compressive_strength
    C  = rheology.ice_compaction_hardening

    E₀ = rheology.undamaged_young_modulus
    η₀ = rheology.undamaged_bulk_viscosity
    α  = rheology.bulk_viscosity_exponent
    
    launch!(architecture(model.grid), model.grid, :xy, 
            _compute_initial_parameters!, P, E, η, P★, C, E₀, η₀, α, h, ℵ)

    fill_halo_regions!(P, model.clock, fields(model))
    
    return nothing
end

# The parameterization for an `ExplicitViscoPlasticRheology`
@kernel function _compute_initial_parameters!(P, E, η, P★, C, E₀, η₀, α, h, ℵ)
    i, j = @index(Global, NTuple)
    
    @inbounds begin 
        K₁ = exp(- C * (1 - ℵ[i, j, 1]))
        K₂ = exp(- α * C * (1 - ℵ[i, j, 1]))
        P[i, j, 1] = P★ * h[i, j, 1] * K₁
        E[i, j, 1] = E₀ * K₁ 
        η[i, j, 1] = η₀ * K₂
    end
end

# Specific compute stresses for the BBM rheology
function compute_stresses!(model, solver, rheology::BrittleBinghamMaxwellRheology, Δt) 
    
    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ℵ  = model.ice_concentration
    ρᵢ = model.ice_density

    substeps = solver.substeps
    stepping_coefficient = solver.substepping_coefficient
    
    u, v = model.velocities
    launch!(arch, grid, :xyz, _compute_bbm_stresses!, rheology, grid, u, v, h, ℵ, ρᵢ, Δt, substeps, stepping_coefficient)

    return nothing
end

@kernel function _compute_bbm_stresses!(rheology, grid, u, v, h, ℵ, ρᵢ, Δt, substeps, stepping_coefficient)
    i, j = @index(Global, NTuple)

    d   = rheology.damage
    P   = rheology.ice_strength
    E   = rheology.young_modulus
    η   = rheology.bulk_viscosity
    σ₁₁ = rheology.σ₁₁
    σ₂₂ = rheology.σ₂₂
    σ₁₂ = rheology.σ₁₂

    @inbounds E₀ = E[i, j, 1]
    @inbounds η₀ = η[i, j, 1]

    λ = η₀ / E₀

    den = 1 / (λ + Δt * (1 - P̃))

    @inbounds σ₁₁′ = λ * (Δt * E * stiffness + σ₁₁[i, j, 1]) * den
    @inbounds σ₂₂′ = λ * (Δt * E * stiffness + σ₂₂[i, j, 1]) * den
    @inbounds σ₁₂′ = λ * (Δt * E * stiffness + σ₁₂[i, j, 1]) * den

    # TODO: finish the rheology
end