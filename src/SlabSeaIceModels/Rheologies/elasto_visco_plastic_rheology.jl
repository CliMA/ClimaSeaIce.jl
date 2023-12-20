#=
struct ElastoViscoPlasticRheology{X, Y, XY} <: AbstractRheology
    substeps :: Int
    Pₚ :: P
    σ₁₁ :: X
    σ₂₂ :: Y
    σ₁₂ :: XY
end

# This will be the implicit step whenever rheology is implemented!
function step_momentum!(model, rheology::ElastoViscoPlastic, Δt, χ)

    u, v = model.velocities
    
    # Alternating leap-frog for Coriolis
    for step in 1:rheology.substeps
        solve_for_stresses!(model, rheology, Δt)
        if iseven(model.clock.iteration)
            solve_for_u!(u, model, rheology, Δt)
            solve_for_v!(v, model, rheology, Δt)
        else
            solve_for_v!(v, model, rheology, Δt)
            solve_for_u!(u, model, rheology, Δt)
        end
    end

    return nothing
end
=#
