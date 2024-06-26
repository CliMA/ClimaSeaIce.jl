# Free drift!

compute_stresses!(model, ::Nothing, args...) = nothing
initialize_rheology!(model, ::Nothing) = nothing

@inline x_internal_stress_divergence(i, j, grid, ::Nothing) = zero(grid)
@inline y_internal_stress_divergence(i, j, grid, ::Nothing) = zero(grid)