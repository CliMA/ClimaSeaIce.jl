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

# Internal stress divergence terms in the velocity tendencies
@inline x_internal_stress_divergenceᶠᶜᶜ(i, j, k, grid, ::Nothing, fields) = zero(grid)
@inline y_internal_stress_divergenceᶠᶜᶜ(i, j, k, grid, ::Nothing, fields) = zero(grid)
@inline x_internal_stress_divergenceᶜᶠᶜ(i, j, k, grid, ::Nothing, fields) = zero(grid)
@inline y_internal_stress_divergenceᶜᶠᶜ(i, j, k, grid, ::Nothing, fields) = zero(grid)

# Additional tendency terms specific to a certain rheology
@inline rheology_specific_numerical_terms_xᶠᶜᶜ(i, j, k, grid, ::Nothing, args...) = zero(grid) 
@inline rheology_specific_numerical_terms_yᶠᶜᶜ(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline rheology_specific_numerical_terms_xᶜᶠᶜ(i, j, k, grid, ::Nothing, args...) = zero(grid) 
@inline rheology_specific_numerical_terms_yᶜᶠᶜ(i, j, k, grid, ::Nothing, args...) = zero(grid)

# Filling the halo regions of the stress components (held in `fields`)
fill_stresses_halo_regions!(fields, dgrid, ::Nothing, args...) = nothing

# Computation of the stresses within the substepping (before the momentum substep)
compute_stresses!(model, solver, ::Nothing, args...) = nothing
