using Oceananigans.Utils: Time
using Oceananigans.Fields: flattened_unique_values
using Oceananigans.OutputReaders: extract_field_time_series, update_field_time_series!

import Oceananigans.Models: update_model_field_time_series!

function step_tracers!(model::SIM, Δt, substep)
    grid = model.grid
    arch = architecture(grid)

    h    = model.ice_thickness
    ℵ    = model.ice_concentration
    hmin = model.ice_consolidation_thickness
    tracers = model.tracers

    Gⁿ = model.timestepper.Gⁿ
    G⁻ = model.timestepper.G⁻

    α, β = timestepping_coefficients(model.timestepper, substep)
    
    launch!(arch, grid, :xyz, _step_tracers!, h, ℵ, hmin, tracers, Gⁿ, G⁻, Δt, α, β)

    return nothing
end

# Thickness and concentration are updated
# We compute hⁿ⁺¹ and ℵⁿ⁺¹ in the same kernel to account for ridging: 
# if ℵ > 1, we reset the concentration to 1 and adjust the thickness 
# to conserve the total ice volume in the cell.
@kernel function _step_tracers!(h, ℵ, hmin, tracers, Gⁿ, G⁻, Δt, α, β)
    i, j, k = @index(Global, NTuple)

    Ghⁿ = Gⁿ.h
    Gℵⁿ = Gⁿ.ℵ
    
    Gh⁻ = G⁻.h
    Gℵ⁻ = G⁻.ℵ

    # Update ice thickness, clipping negative values
    @inbounds begin
        h⁻ = hmin[i, j, k]
        h⁺ = h[i, j, k] + Δt * (α * Ghⁿ[i, j, k] + β * Gh⁻[i, j, k])
        h⁺ = max(zero(h⁺), h⁺)
        
        ice_covered = h⁺ > 0
        ℵ⁺ = ℵ[i, j, k] + Δt * (α * Gℵⁿ[i, j, k] + β * Gℵ⁻[i, j, k])
        ℵ⁺ = max(zero(ℵ⁺), ℵ⁺) * ice_covered

        # If h < hmin we reset the thickness to zero and adjust the concentration accordingly
        # to maintain a constant ice volume
        ht = ifelse(h⁺ < h⁻, h⁻, h⁺)
        ℵt = ifelse(ht < h⁻, ℵt * (h⁻ - h⁺) / h⁺, ℵt)
        ℵt = ifelse(ht == 0, zero(ℵt), ℵt)

        # If ℵ > 1, we reset the concentration to 1 and increase the thickness accordingly
        # to maintain a constant ice volume
        ht = ifelse(ℵt > 1, ht * ℵt, ht)
        ℵt = ifelse(ℵt > 1, one(ℵt), ℵt)

        ℵ[i, j, k] = ℵt
        h[i, j, k] = ht
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

    update_model_field_time_series!(model, model.clock)

    return nothing
end

function update_model_field_time_series!(model::SeaIceModel, clock::Clock)
    time = Time(clock.time)

    possible_fts = (model.tracers, model.external_heat_fluxes, model.external_momentum_stresses)
    time_series_tuple = extract_field_time_series(possible_fts)
    time_series_tuple = flattened_unique_values(time_series_tuple)

    for fts in time_series_tuple
        update_field_time_series!(fts, time)
    end

    return nothing
end
