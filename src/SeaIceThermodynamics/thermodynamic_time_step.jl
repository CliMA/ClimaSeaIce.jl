using Oceananigans.Architectures: architecture
using Oceananigans.Utils
using KernelAbstractions: @kernel, @index

thermodynamic_time_step!(model, ::Nothing, Δt) = nothing

function thermodynamic_time_step!(model, ::SlabSeaIceThermodynamics, Δt)
    grid = model.grid
    arch = architecture(grid)
    
    launch!(arch, grid, :xy,
            _slab_thermodynamic_time_step!,
            model.ice_thickness,
            model.ice_concentration,
            grid, Δt,
            model.clock,
            model.ice_consolidation_thickness,
            model.ice_thermodynamics,
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
@kernel function _slab_thermodynamic_time_step!(ice_thickness,
                                           ice_concentration,
                                           grid,
                                           Δt,
                                           clock,
                                           ice_consolidation_thickness,
                                           ice_thermodynamics,
                                           top_external_heat_flux,
                                           bottom_external_heat_flux,
                                           model_fields)

    i, j = @index(Global, NTuple)
    
    @inbounds hⁿ = ice_thickness[i, j, 1]
    @inbounds ℵⁿ = ice_concentration[i, j, 1]
    @inbounds hᶜ = ice_consolidation_thickness[i, j, 1]

    # Total volume tendency
    ∂t_V = thermodynamic_tendency(i, j, 1, grid,
                                  ice_thermodynamics,
                                  ice_thickness,
                                  ice_concentration,
                                  ice_consolidation_thickness,
                                  top_external_heat_flux,
                                  bottom_external_heat_flux,
                                  clock, model_fields)

    # ice volume at timestep n+1
    Vⁿ⁺¹ = hⁿ * ℵⁿ + Δt * ∂t_V

    # Adjust the ice volume to zero
    Vⁿ⁺¹ = max(zero(Vⁿ⁺¹), Vⁿ⁺¹)

    # We recalculate the actual volume derivative, after accounting for the
    # volume adjustment (the ice cannot produce more melt than its actual volume!)
    ∂t_V = (Vⁿ⁺¹ - hⁿ * ℵⁿ) / Δt

    # We parameterize the evolution of ice thickness and concentration
    # (i.e. lateral vs vertical growth) following Hibler (1979)
    ∂t_V_freezing = max(∂t_V, zero(ℵⁿ))
    ∂t_V_melting  = min(∂t_V, zero(ℵⁿ))

    ∂t_ℵᶠ = (1 - ℵⁿ) /  hᶜ * ∂t_V_freezing
    ∂t_ℵᵐ =      ℵⁿ  / 2hⁿ * ∂t_V_melting
    
    # Update ice thickness and concentration
    ℵ⁺ = ℵⁿ + Δt * (∂t_ℵᶠ + ∂t_ℵᵐ)
    h⁺ = Vⁿ⁺¹ / ℵ⁺

    # Treat pathological cases
    ℵ⁺ = max(zero(ℵ⁺), ℵ⁺)
    h⁺ = ifelse(ℵ⁺ ≤ 0, zero(h⁺), h⁺)

    # No volume change
    ℵ⁺ = ifelse(∂t_V == 0, ℵⁿ, ℵ⁺)
    h⁺ = ifelse(∂t_V == 0, hⁿ, h⁺)
    ℵ⁺ = ifelse(h⁺ == 0, zero(ℵ⁺), ℵ⁺) # reset the concentration if there is no sea-ice

    # Ridging caused by the thermodynamic step
    @inbounds ice_concentration[i, j, 1] = ifelse(ℵ⁺ > 1, one(ℵ⁺), ℵ⁺)
    @inbounds ice_thickness[i, j, 1]     = ifelse(ℵ⁺ > 1,  h⁺ * ℵ⁺, h⁺)
end
