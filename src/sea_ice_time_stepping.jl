function step_tracers!(model::SIM, Δt, substep)
    grid = model.grid
    arch = architecture(grid)

    h  = model.ice_thickness
    ℵ  = model.ice_concentration
    tracers = model.tracers

    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻

    α, β = timestepping_coefficients(model.timestepper, substep)
    launch!(arch, grid, :xyz, _step_tracers!, h, ℵ, tracers, Gⁿ, G⁻, Δt, α, β)

    return nothing
end

# Thickness and concentration are updated
# We compute hⁿ⁺¹ and ℵⁿ⁺¹ in the same kernel to account for ridging: 
# if ℵ > 1, we reset the concentration to 1 and adjust the thickness 
# to conserve the total ice volume in the cell.
@kernel function _step_tracers!(h, ℵ, tracers, Gⁿ, G⁻, Δt, α, β)
    i, j, k = @index(Global, NTuple)

    Ghⁿ = Gⁿ.h
    Gℵⁿ = Gⁿ.ℵ
    Gh⁻ = G⁻.h
    Gℵ⁻ = G⁻.ℵ

    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁺ = h[i, j, k] + Δt * (α * Ghⁿ[i, j, k] + β * Gh⁻[i, j, k])
        h⁺ = max(0, h⁺)

        ℵ⁺ = ℵ[i, j, k] + Δt * (α * Gℵⁿ[i, j, k] + β * Gℵ⁻[i, j, k])
        ℵ⁺ = max(0, ℵ⁺)
        
        # Ridging! if ℵ > 1, we reset the concentration to 1 and increase the thickness accordingly
        # to maintain a constant ice volume
        h[i, j, k] = ifelse(ℵ⁺ > 1, h⁺ * ℵ⁺, h⁺)
        ℵ[i, j, k] = ifelse(ℵ⁺ > 1, one(ℵ⁺), ℵ⁺)
    end 
end

function store_tendencies!(model::SIM) 

    grid = model.grid
    arch = architecture(grid)
    Nx, Ny, _ = size(grid)

    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻
    Nt = length(Gⁿ)

    params = KernelParameters((Nx, Ny, Nt), (0, 0, 0))
    launch!(arch, model.grid, params, _store_tendencies!, G⁻, Gⁿ)

    return nothing
end

@kernel function _store_tendencies!(G⁻, Gⁿ) 
    i, j, n = @index(Global, NTuple)
    @inbounds G⁻[n][i, j, 1] = Gⁿ[n][i, j, 1]
end

const RK3SeaIceModel = SeaIceModel{<:Any, <:Any, <:Any, <:RungeKutta3TimeStepper}

function time_step!(model::RK3SeaIceModel, Δt; callbacks = [])

    # Be paranoid and update state at iteration 0, in case run! is not used:
    model.clock.iteration == 0 && update_state!(model)

    γ¹ = model.timestepper.γ¹
    γ² = model.timestepper.γ²
    γ³ = model.timestepper.γ³

    ζ² = model.timestepper.ζ²
    ζ³ = model.timestepper.ζ³

    first_stage_Δt  = γ¹ * Δt
    second_stage_Δt = (γ² + ζ²) * Δt
    third_stage_Δt  = (γ³ + ζ³) * Δt

    #
    # First stage
    #

    compute_tendencies!(model)
    step_tracers!(model, Δt, 1)
    step_momentum!(model, model.ice_dynamics, Δt)
    store_tendencies!(model)

    tick!(model.clock, first_stage_Δt)
    update_state!(model)

    #
    # Second stage
    #

    compute_tracer_tendencies!(model)
    step_tracers!(model, Δt, 2)
    step_momentum!(model, model.ice_dynamics, Δt)
    store_tendencies!(model)

    tick!(model.clock, second_stage_Δt)
    update_state!(model)

    #
    # Third stage
    #

    compute_tracer_tendencies!(model)
    step_tracers!(model, Δt, 3)
    step_momentum!(model, model.ice_dynamics, Δt)
    store_tendencies!(model)

    tick!(model.clock, third_stage_Δt)
    update_state!(model)

    return nothing
end

function timestepping_coefficients(ts::RungeKutta3TimeStepper, substep) 
    if substep == 1
        return ts.γ¹, zero(ts.γ¹)
    elseif substep == 2
        return ts.γ², ts.ζ²
    elseif substep == 3
        return ts.γ³, ts.ζ³
    end
end

function update_state!(model::SIM)
    foreach(prognostic_fields(model)) do field
        mask_immersed_field!(field)
    end

    fill_halo_regions!(prognostic_fields(model), model.clock, fields(model))

    return nothing
end


