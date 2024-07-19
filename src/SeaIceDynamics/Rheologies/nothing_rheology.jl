####
#### Free drift!
####

# This file shows all the functions that need to be extended 
# when defining a new rheology compatible with the `ExplicitMomentumSolver`

# Fields needed to compute stresses (and stresses themselves)
required_auxiliary_fields(::Nothing) = nothing

# All initializations required before the substepping begins
initialize_rheology!(model, ::Nothing) = nothing

# Computation of the stresses within the substepping (before the momentum substep)
compute_stresses!(model, solver, ::Nothing, args...) = nothing

# Internal stress divergence terms in the velocity tendencies
@inline ∂ⱼ_σ₁ⱼ(i, j, k, grid, ::Nothing) = zero(grid)
@inline ∂ⱼ_σ₂ⱼ(i, j, k, grid, ::Nothing) = zero(grid)

# Additional tendency terms specific to a certain rheology
@inline rheology_specific_forcing_x(i, j, k, grid, args...) = zero(grid) 
@inline rheology_specific_forcing_y(i, j, k, grid, args...) = zero(grid)
