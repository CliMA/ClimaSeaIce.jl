
struct BrittleBinghamMaxwellRheology{S1, S2, S3, P, FT} <: AbstractRheology
    σ₁₁   :: S1 # internal stress xx
    σ₂₂   :: S2 # internal stress yy
    σ₁₂   :: S3 # internal stress xy
    ice_strength :: P # field containing the precomputed ice strength
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
end

function BrittleBinghamMaxwellRheology(grid::AbstractGrid; 
                                       ice_compressive_strength = 27500, 
                                       ice_compaction_hardening = 20)
    σ₁₁ = CenterField(grid)
    σ₂₂ = CenterField(grid)
    σ₁₂ = Field{Face, Face, Center}(grid)

    P  = CenterField(grid)
    FT = eltype(grid)
    return BrittleBinghamMaxwellRheology(σ₁₁, σ₂₂, σ₁₂, P, 
                                         convert(FT, ice_compressive_strength), 
                                         convert(FT, ice_compaction_hardening))
end

# Extend the `adapt_structure` function for the BrittleBinghamMaxwellRheology
Adapt.adapt_structure(to, r::BrittleBinghamMaxwellRheology) = 
    BrittleBinghamMaxwellRheology(Adapt.adapt(to, r.σ₁₁),
                                  Adapt.adapt(to, r.σ₂₂),
                                  Adapt.adapt(to, r.σ₁₂),
                                  Adapt.adapt(to, r.ice_strength),
                                  Adapt.adapt(to, r.ice_compressive_strength),
                                  Adapt.adapt(to, r.ice_compaction_hardening))

# Specific compute stresses for the EVP rheology
function compute_stresses!(model, solver, rheology::ExplicitViscoPlasticRheology, Δt) 

end

