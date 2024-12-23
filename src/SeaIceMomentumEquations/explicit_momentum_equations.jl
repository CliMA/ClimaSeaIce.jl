using ClimaSeaIce: timestepping_coefficients

struct ExplicitMomentumEquation{C, R, A}
    coriolis :: C
    rheology :: R
    auxiliary_fields :: A
end

function ExplicitMomentumEquation(grid; 
                                  coriolis = nothing,
                                  rheology = nothing,
                                  auxiliary_fields = NamedTuple())

    auxiliary_fields = merge(auxiliary_fields, required_auxiliary_fields(rheology, grid))

    return ExplicitMomentumEquation(coriolis, rheology, auxiliary_fields)
end

# Simple explicit stepping of the momentum equations
function step_momentum!(model, ::ExplicitMomentumEquation, Δt, stage)
    grid = model.grid
    arch = architecture(grid)

    u, v = model.velocities
    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻

    α, β = timestepping_coefficients(model.timestepper, stage)
    
    launch!(arch, grid, :xyz, _step_velocities!, u, v, Gⁿ, G⁻, Δt, α, β)

    return nothing
end

@kernel function _step_velocities!(u, v, Gⁿ, G⁻, Δt, α, β)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        u[i, j, k] += Δt * (α * Gⁿ.u[i, j, k] + β * G⁻.u[i, j, k])
        v[i, j, k] += Δt * (α * Gⁿ.v[i, j, k] + β * G⁻.v[i, j, k])
    end 
end