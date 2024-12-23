module Rheologies

export ViscousRheology

using Oceananigans

# The only functions we need to extend
import ClimaSeaIce: ∂ⱼ_σ₁ⱼ, ∂ⱼ_σ₂ⱼ, required_auxiliary_fields

# Nothing rheology
required_auxiliary_fields(rheology, grid) = NamedTuple()
@inline ∂ⱼ_σ₁ⱼ(i, j, k, grid, rheology, fields) = zero(grid)
@inline ∂ⱼ_σ₂ⱼ(i, j, k, grid, rheology, fields) = zero(grid)

include("viscous_rheology.jl")

end