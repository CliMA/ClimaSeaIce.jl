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

    FT = eltype(χ)
    Δt = convert(FT, Δt)

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

function update_state!(model::SIM, callbacks = nothing)
    
    foreach(prognostic_fields(model)) do field
        mask_immersed_field!(field)
    end

    fill_halo_regions!(prognostic_fields(model), model.clock, fields(model))

    return nothing
end
