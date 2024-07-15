####
#### Free drift!
####

# This file shows all the functions that need to be extended 
# when defining a new rheology compatible with the `ExplicitMomentumSolver`

required_auxiliary_fields(::Nothing) = nothing

compute_stresses!(model, solver, ::Nothing, args...) = nothing
initialize_rheology!(model, ::Nothing) = nothing

@inline x_internal_stress_divergence(i, j, k, grid, ::Nothing) = zero(grid)
@inline y_internal_stress_divergence(i, j, k, grid, ::Nothing) = zero(grid)

@inline rheology_specific_numerical_terms_x(i, j, k, grid, args...) = zero(grid) 
@inline rheology_specific_numerical_terms_y(i, j, k, grid, args...) = zero(grid)

fill_stresses_halo_regions!(fields, ::Nothing, args...) = nothing