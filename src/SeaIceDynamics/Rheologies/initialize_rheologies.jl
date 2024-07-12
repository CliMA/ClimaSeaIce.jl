
# Extend `fill_halo_regions!` for the AbstractRheology
function fill_halo_regions!(rheology::AbstractRheology, args...)
    fill_halo_regions!(rheology.σ₁₁, args...)
    fill_halo_regions!(rheology.σ₂₂, args...)
    fill_halo_regions!(rheology.σ₁₂, args...)

    return nothing
end

# Extend `mask_immersed_field!` for the AbstractRheology
function mask_immersed_field!(rheology::AbstractRheology)
    mask_immersed_field!(rheology.σ₁₁)
    mask_immersed_field!(rheology.σ₂₂)
    mask_immersed_field!(rheology.σ₁₂)

    return nothing
end

"""
    initialize_rheology!(model, rheology::AbstractRheology)

Initialize the elasto-visco-plastic rheology.
In this step we calculate the ice strength given the ice mass (thickness and concentration).
"""
function initialize_rheology!(model, rheology::AbstractRheology)
    h = model.ice_thickness
    ℵ = model.concentration

    P  = rheology.ice_strength
    P★ = rheology.ice_compressive_strength
    C  = rheology.ice_compaction_hardening
    
    launch!(architecture(model.grid), model.grid, :xy, _compute_ice_strength!, P, P★, C, h, ℵ)

    fill_halo_regions!(P, model.clock, fields(model))
    
    return nothing
end

# The parameterization for an `ExplicitViscoPlasticRheology`
@kernel function _compute_ice_strength!(P, P★, C, h, ℵ)
    i, j = @index(Global, NTuple)
    @inbounds P[i, j, 1] = @inbounds P★ * h[i, j, 1] * exp(- C * (1 - ℵ[i, j, 1])) 
end