
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
