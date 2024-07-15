# Free drift!
required_auxiliary_fields(grid, ::Nothing, ::CgridDynamics) = nothing
required_auxiliary_fields(grid, ::Nothing, ::EgridDynamics) = (û = YFaceField(grid), v̂ = XFaceField(grid))

compute_stresses!(model, solver, ::Nothing, args...) = nothing
initialize_rheology!(model, ::Nothing) = nothing

@inline x_internal_stress_divergence(i, j, k, grid, ::Nothing) = zero(grid)
@inline y_internal_stress_divergence(i, j, k, grid, ::Nothing) = zero(grid)

@inline rheology_specific_numerical_terms_x(i, j, k, grid, args...) = zero(grid) 
@inline rheology_specific_numerical_terms_y(i, j, k, grid, args...) = zero(grid)