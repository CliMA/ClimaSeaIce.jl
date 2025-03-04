using Oceananigans.Architectures: architecture
using Oceananigans.Utils
using KernelAbstractions: @kernel, @index

thermodynamic_step!(model, ::Nothing, Δt) = nothing

function thermodynamic_step!(model, ::SlabSeaIceThermodynamics, Δt)
    grid = model.grid
    arch = architecture(grid)
    
    launch!(arch, grid, :xy,
            _slab_thermodynamic_step!,
            model.ice_thickness,
            model.ice_concentration,
            grid, Δt,
            model.clock,
            model.ice_consolidation_thickness,
            model.thermodynamics,
            model.external_heat_fluxes.top,
            model.external_heat_fluxes.bottom,
            fields(model))


    return nothing
end

# The thermodynamic step is computed in a single kernel following:
#
# ∂t_V = ∂t_h * ℵ + h * ∂t_ℵ
#
# Therefore:
#    
#                     h⁺ * ℵ⁺ - hⁿ * ℵⁿ  
#  ∂t_V = Gᴸ + Gⱽ  = ------------------- 
#                             Δt               
#      
# The two will be adjusted conservatively after the thermodynamic step to ensure that ℵ ≤ 1.
@kernel function _slab_thermodynamic_step!(ice_thickness,
                                           ice_concentration,
                                           grid,
                                           Δt,
                                           clock,
                                           ice_consolidation_thickness,
                                           thermodynamics,
                                           top_external_heat_flux,
                                           bottom_external_heat_flux,
                                           model_fields)

    i, j = @index(Global, NTuple)
    
    @inbounds hⁿ = ice_thickness[i, j, 1]
    @inbounds ℵⁿ = ice_concentration[i, j, 1]

    Gⱽ = vertical_growth(i, j, 1, grid,
                         thermodynamics,
                         ice_thickness,
                         ice_concentration,
                         ice_consolidation_thickness,
                         top_external_heat_flux,
                         bottom_external_heat_flux,
                         clock, model_fields)

    Gᴸ = lateral_growth(i, j, 1, grid, thermodynamics, bottom_external_heat_flux, clock, model_fields)

    # Total volume tendency
    ∂t_V = Gᴸ + Gⱽ

    # ice volume at timestep n+1
    Vⁿ⁺¹ = hⁿ * ℵⁿ + Δt * ∂t_V

    # Adjust the ice volume to zero
    Vⁿ⁺¹ = max(zero(Vⁿ⁺¹), Vⁿ⁺¹)

    # If Vⁿ⁺¹ == 0 the ice has melted completely, and we set both hⁿ⁺¹ and ℵⁿ⁺¹ to zero
    # Otherwise, if the volume is positive, we adjust the thickness and concentration conservatively
    # To account for this, we recalculate the actual volume derivative
    ∂t_V = (Vⁿ⁺¹ - hⁿ * ℵⁿ) / Δt

    # Simple explicit step, we assume lateral growth 
    # (at the beginning) contributes only to the ice concentration
    ℵ⁺ = ℵⁿ + Δt * Gᴸ / hⁿ * (hⁿ > 0)
    ℵ⁺ = max(zero(ℵ⁺), ℵ⁺) # Concentration cannot be negative, clip it up

    # The concentration derivative
    ∂t_ℵ = (ℵ⁺ - ℵⁿ) / Δt

    # Adjust the thickness accordingly
    h⁺ = hⁿ + Δt * (∂t_V - hⁿ * ∂t_ℵ) / ℵ⁺ * (ℵ⁺ > 0)

    # Ridging and rafting caused by the thermodynamic step
    h⁺ = max(zero(h⁺), h⁺) # Thickness cannot be negative, clip it up
    ℵ⁺ = ifelse(Vⁿ⁺¹ == 0, zero(ℵ⁺), ℵ⁺)
    h⁺ = ifelse(Vⁿ⁺¹ == 0, zero(h⁺), h⁺)
    
    @inbounds ice_concentration[i, j, 1] = ifelse(ℵ⁺ > 1, one(ℵ⁺), ℵ⁺)
    @inbounds ice_thickness[i, j, 1]     = ifelse(ℵ⁺ > 1,  h⁺ * ℵ⁺, h⁺)
end