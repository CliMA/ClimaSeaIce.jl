using Oceananigans.Utils
using ClimaSeaIce: timestepping_coefficients

const ExplicitMomentumEquation = SeaIceMomentumEquation{<:ExplicitSolver}

# Simple explicit stepping of the momentum equations
function step_momentum!(model, ::ExplicitMomentumEquation, Δt, stage)
    grid = model.grid
    arch = architecture(grid)

    u, v = model.velocities
    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻

    initialize_rheology!(model, ice_dynamics.rheology)
    compute_stresses!(model, ice_dynamics, ice_dynamics.rheology, Δt)

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

# Compute the tendencies for the explicit momentum equations
function compute_momentum_tendencies!(model, ::ExplicitMomentumEquation)
    
    ice_dynamics = model.ice_dynamics
    grid = model.grid

    args = (model.clock, 
            model.velocities, 
            ice_dynamics.coriolis, 
            ice_dynamics.rheology, 
            ice_dynamics.auxiliary_fields, 
            model.ice_thickness, 
            model.ice_concentration, 
            model.ice_density)

    u_top_stress = model.external_momentum_stresses.top.u
    v_top_stress = model.external_momentum_stresses.top.v

    u_bottom_stress = model.external_momentum_stresses.bottom.u
    v_bottom_stress = model.external_momentum_stresses.bottom.v

    Gu = model.timestepper.Gⁿ.u
    Gv = model.timestepper.Gⁿ.v

    launch!(architecture(grid), grid, :xy, _compute_velocity_tendencies!, Gu, Gv, grid, args,
            u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress)

    return nothing
end

@kernel function _compute_velocity_tendencies!(Gu, Gv, grid, args, u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress)
    i, j = @index(Global, NTuple)
    @inbounds Gu[i, j, 1] = u_velocity_tendency(i, j, grid, args..., u_top_stress, u_bottom_stress)
    @inbounds Gv[i, j, 1] = v_velocity_tendency(i, j, grid, args..., v_top_stress, v_bottom_stress)
end