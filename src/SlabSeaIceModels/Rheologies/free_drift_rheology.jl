"""
    FreeDriftRheology

a Type for free drift rheologies, where ``∇ ⋅ σ = 0 ``
"""
struct FreeDriftRheology <: AbstractExplicitRheology
    substeps :: Int
end

compute_stresses!(model, ::FreeDriftRheology, Δτ) = nothing

@inline x_stress_divergence(i, j, grid, ::FreeDriftRheology) = zero(grid)
@inline y_stress_divergence(i, j, grid, ::FreeDriftRheology) = zero(grid)