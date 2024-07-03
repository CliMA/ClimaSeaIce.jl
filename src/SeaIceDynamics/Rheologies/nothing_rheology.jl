# Free drift!

compute_stresses!(model, solver, ::Nothing, args...) = nothing
initialize_rheology!(model, ::Nothing) = nothing

@inline x_internal_stress_divergence(i, j, k, grid, ::Nothing) = zero(grid)
@inline y_internal_stress_divergence(i, j, k, grid, ::Nothing) = zero(grid)