####
#### Nothing rheology a.k.a. "Free drift" (fallback implementation)
####

# This file shows all the functions that need to be extended 
# when defining a new rheology compatible with the `ExplicitDynamics`

# Fields needed to compute stresses (and stresses themselves)
required_auxiliary_fields(grid, rheology) = NamedTuple()

# All initializations required before the substepping begins
initialize_rheology!(model, rheology) = nothing

# Computation of the stresses within the substepping (before the momentum substep)
compute_stresses!(model, solver, rheology, args...) = nothing

# Internal stress divergence terms in the velocity tendencies
@inline ∂ⱼ_σ₁ⱼ(i, j, k, grid, rheology, fields) = zero(grid)
@inline ∂ⱼ_σ₂ⱼ(i, j, k, grid, rheology, fields) = zero(grid)

# Additional tendency terms specific to a certain rheology
@inline rheology_specific_forcing_x(i, j, k, grid, args...) = zero(grid) 
@inline rheology_specific_forcing_y(i, j, k, grid, args...) = zero(grid)

# Fallback for no specific scheme => α = β = substeps 
@inline rheology_substeps(i, j, k, grid, rheology, substeps, fields) = substeps 
