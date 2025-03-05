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
    @inbounds hᶜ = ice_consolidation_thickness[i, j, 1]

    # Total volume tendency
    ∂t_V = thermodynamic_tendency(i, j, 1, grid,
                                  thermodynamics,
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

    # There are 6 typical cases depending on the sign of ∂t_V, ℵⁿ and hⁿ
    #
    # - ∂t_V ≥ 0 : the ice is growing * easy case to handle *
    # ├── ℵⁿ == 0 -> new ice is forming (we only increase ℵⁿ, post increase ridging will adjust h)
    # └── ℵⁿ > 0  -> ice is growing (we increase both ℵ and h based on ℵⁿ)
    #
    # - ∂t_V < 0 : the ice is melting * tricky case to handle *
    # ├── ℵⁿ == 0 -> not possible! (there is no ice to begin with)
    # ├── h > hᶜ  -> ice is melting (we decrease both ℵ and h based on ℵⁿ)
    # └── h ≤ hᶜ  -> unconsolidated ice is melting (we only decrease ℵⁿ)

    freezing     = ∂t_V > 0  
    consolidated = hⁿ > hᶜ
    open_ocean   = ℵⁿ == 0

    # Freezing and melting cases:
    hᶠ = hⁿ + Δt * ∂t_V * ℵⁿ * !open_ocean
    hᵐ = hⁿ + Δt * ∂t_V * ℵⁿ * consolidated
    hᵐ = max(hᵐ, zero(hᵐ))
    h⁺ = ifelse(freezing, hᶠ, hᵐ)

    # There is also a very particular case when the ice is nucleating corresponding
    # h⁺ == 0 and freezing. In this case, we set h = hᶜ and advance ℵ accordingly
    h⁺ = ifelse(freezing & (h⁺ == 0), hᶜ, h⁺)
    ℵ⁺ = Vⁿ⁺¹ / h⁺
    
    # No volume change
    ℵ⁺ = ifelse(∂t_V == 0, ℵⁿ, ℵ⁺)
    h⁺ = ifelse(∂t_V == 0, hⁿ, h⁺)

    # Ridging caused by the thermodynamic step
    @inbounds ice_concentration[i, j, 1] = ifelse(ℵ⁺ > 1, one(ℵ⁺), ℵ⁺)
    @inbounds ice_thickness[i, j, 1]     = ifelse(ℵ⁺ > 1,  h⁺ * ℵ⁺, h⁺)
end