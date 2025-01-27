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
        h⁺ = max(zero(h⁺), h⁺)

        ice_covered = h⁺ > 0
        ℵ⁺ = ℵ[i, j, k] + Δt * (α * Gℵⁿ[i, j, k] + β * Gℵ⁻[i, j, k])
        ℵ⁺ = max(zero(ℵ⁺), ℵ⁺) * ice_covered
        
        # Ridging! if ℵ > 1, we reset the concentration to 1 and increase the thickness accordingly
        # to maintain a constant ice volume
        h[i, j, k] = ifelse(ℵ⁺ > 1, h⁺ * ℵ⁺, h⁺)
        h[i, j, k] = ifelse(ℵ⁺ == 0, zero(h⁺), h⁺)
        ℵ[i, j, k] = ifelse(ℵ⁺ > 1, one(ℵ⁺), ℵ⁺)
    end 
end


function store_tendencies!(model::SIM) 

    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻
    Nt = length(Gⁿ)

    for n in 1:Nt
        parent(G⁻[n]) .= parent(Gⁿ[n])
    end

    return nothing
end

function update_state!(model::SIM)
    
    foreach(prognostic_fields(model)) do field
        mask_immersed_field!(field)
    end

    fill_halo_regions!(prognostic_fields(model), model.clock, fields(model))

    return nothing
end

