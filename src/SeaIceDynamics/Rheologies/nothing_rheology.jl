####
#### Free drift!
####

# This file shows all the functions that an "explicit" rheology should 
# extend to work correctly with an `ExplicitMomentumSolver`

# Fields needed to compute stresses (and stresses themselves)
required_auxiliary_fields(grid, ::Nothing, ::CGridDynamics) = nothing
required_auxiliary_fields(grid, ::Nothing, ::EGridDynamics) = (û = YFaceField(grid), v̂ = XFaceField(grid))

# All initializations required before the substepping begins
initialize_rheology!(model, ::Nothing) = nothing

# Computation of the stresses within the substepping (before the momentum substep)
compute_stresses!(model, solver, ::Nothing, args...) = nothing

# Internal stress divergence terms in the velocity tendencies
@inline ∂ⱼ_σ₁ⱼᶠᶜᶜ(i, j, k, grid, ::Nothing) = zero(grid)
@inline ∂ⱼ_σ₂ⱼᶠᶜᶜ(i, j, k, grid, ::Nothing) = zero(grid)
@inline ∂ⱼ_σ₁ⱼᶜᶠᶜ(i, j, k, grid, ::Nothing) = zero(grid)
@inline ∂ⱼ_σ₂ⱼᶜᶠᶜ(i, j, k, grid, ::Nothing) = zero(grid)

# Additional tendency terms specific to a certain rheology
@inline rheology_specific_forcing_xᶠᶜᶜ(i, j, k, grid, ::Nothing, args...) = zero(grid) 
@inline rheology_specific_forcing_yᶠᶜᶜ(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline rheology_specific_forcing_xᶜᶠᶜ(i, j, k, grid, ::Nothing, args...) = zero(grid) 
@inline rheology_specific_forcing_yᶜᶠᶜ(i, j, k, grid, ::Nothing, args...) = zero(grid)

