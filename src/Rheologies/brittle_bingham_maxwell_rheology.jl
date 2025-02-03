## The equations are solved in an iterative form following the EVP rheology of
## Kimmritz et al (2016) (https://www.sciencedirect.com/science/article/pii/S1463500317300690)
#
# Where:
# σᵢⱼ(u) = 2η ϵ̇ᵢⱼ + [(ζ - η) * (ϵ̇₁₁ + ϵ̇₂₂) - P / 2] δᵢⱼ
#
struct BrittleBinghamMaxellRheology{FT}
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    yield_curve_eccentricity :: FT # elliptic yield curve eccentricity
    minimum_plastic_stress :: FT # minimum plastic parameter (transitions to viscous behaviour)
    min_substeps :: FT # minimum number of substeps expressed as the dynamic coefficient
    max_substeps :: FT # maximum number of substeps expressed as the dynamic coefficient
end


function BrittleBinghamMaxellRheology(FT::DataType = Float64; 
                                    ice_compressive_strength = 27500, 
                                    ice_compaction_hardening = 20, 
                                    yield_curve_eccentricity = 2, 
                                    minimum_plastic_stress = 2e-9,
                                    min_substeps = 30,
                                    max_substeps = 1000)

    return BrittleBinghamMaxellRheology(convert(FT, ice_compressive_strength), 
                                         convert(FT, ice_compaction_hardening), 
                                         convert(FT, yield_curve_eccentricity),
                                         convert(FT, minimum_plastic_stress),
                                         convert(FT, min_substeps),
                                         convert(FT, max_substeps))
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
    ElastoViscoPlasticRheology(Adapt.adapt(to, r.ice_compressive_strength),
                               Adapt.adapt(to, r.ice_compaction_hardening),
                               Adapt.adapt(to, r.yield_curve_eccentricity),
                               Adapt.adapt(to, r.minimum_plastic_stress),
                               Adapt.adapt(to, r.min_substeps),
                               Adapt.adapt(to, r.max_substeps))

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