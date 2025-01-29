using Oceananigans.Utils
using ClimaSeaIce: timestepping_coefficients

const ExplicitMomentumEquation = SeaIceMomentumEquation{<:ExplicitSolver}
 
#=

This needs to be extended for the new solver

# Compute the tendencies for the explicit momentum equations
function step_momentum!(model, ::ExplicitMomentumEquation, Δt, stage)
    
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

    u, v = model.velocities

    u_forcing = model.forcing.u
    v_forcing = model.forcing.v

    initialize_rheology!(model, ice_dynamics.rheology)
    compute_stresses!(model, ice_dynamics, ice_dynamics.rheology, Δt)

    α, β = timestepping_coefficients(model.timestepper, stage)
    
    launch!(architecture(grid), grid, :xy, _step_velocities!, u, v, Gu, Gv, grid, Δt, α, β, args,
            u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress, u_forcing, v_forcing)

    return nothing
end

@kernel function _compute_velocity_tendencies!(u, v, Gu, Gv, grid, Δt, α, β, args, 
                                               u_top_stress, v_top_stress, u_bottom_stress, v_bottom_stress,
                                               u_forcing, v_forcing)
    i, j = @index(Global, NTuple)
    τuᵢ, Guⁿ = u_velocity_tendency(i, j, grid, args..., u_top_stress, u_bottom_stress, u_forcing)
    τvᵢ, Gvⁿ = v_velocity_tendency(i, j, grid, args..., v_top_stress, v_bottom_stress, v_forcing)

    @inbounds u[i, j, k] = (u[i, j, k] + Δt * (α * Guⁿ + β * Gu[i, j, k])) / (1 + Δt * τuᵢ)
    @inbounds v[i, j, k] = (v[i, j, k] + Δt * (α * Gvⁿ + β * Gv[i, j, k])) / (1 + Δt * τvᵢ)

    # Store directly the tendencies
    @inbounds Gu[i, j, k] = Guⁿ
    @inbounds Gv[i, j, k] = Gvⁿ
end
=#

s