
thermodynamic_step!(model, ::Nothing, Δt) = nothing

function thermodynamic_step!(model, ::SlabSeaIceThermodynamics, Δt)
    grid = model.grid
    arch = architecture(grid)
    
    launch!(arch, grid, :xy,
            _slab_thermodynamic_step!,
            model.ice_thickness,
            model.ice_concentration,
            grid,
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
# We consider the lateral growth as increasing ℵ and the vertical growth as increasing h. 
# This equation is solved semi-implicitly to avoid NaNs coming from the division by quantities that can go to zero.
# Therefore:
#    
#                     hⁿ⁺¹ - hⁿ          ℵⁿ⁺¹ - ℵⁿ       
#  ∂t_V = Gᴸ + Gⱽ  = --------- ⋅ ℵⁿ⁺¹ + --------- ⋅ hⁿ⁺¹
#                       Δt                  Δt           
#      
# Leading to:
#  
#  ℵⁿ⁺¹ - ℵⁿ     Gᴸ    
#  --------- = ----- 
#     Δt        hⁿ⁺¹  
#     
# And
#
#  hⁿ⁺¹ - hⁿ            hⁿ⁺¹ - 1    1
#  --------- = [Gⱽ + Gᴸ --------] -----
#     Δt                 hⁿ⁺¹      ℵⁿ⁺¹  
#
# The two will be adjusted conservatively after the thermodynamic step to ensure that ℵ ≤ 1.
@kernel function _slab_thermodynamic_step!(ice_thickness,
                                           ice_concentration,
                                           grid,
                                           clock,
                                           ice_consolidation_thickness,
                                           thermodynamics,
                                           top_external_heat_flux,
                                           bottom_external_heat_flux,
                                           model_fields)

    i, j = @index(Global, NTuple)
    
    @inbounds hᵢ = ice_thickness[i, j, 1]
    @inbounds ℵᵢ = ice_concentration[i, j, 1]

    Gh⁺ = lateral_growth(i, j, 1, grid,
                           thermodynamics,
                           ice_thickness,
                           ice_concentration,
                           ice_consolidation_thickness,
                           top_external_heat_flux,
                           bottom_external_heat_flux,
                           clock, model_fields)

    GV⁺ = vertical_growth(i, j, 1, grid, thermodynamics, bottom_external_heat_flux, clock, model_fields)

    Gℵ = GV⁺ / hᵢ
    Gh = Gh⁺ * ℵᵢ

    Guh = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, ice_thickness)
    Guℵ = - horizontal_div_Uc(i, j, 1, grid, advection, velocities, ice_concentration)

    @inbounds Gⁿ.h[i, j, 1] = Gh + Guh
    @inbounds Gⁿ.ℵ[i, j, 1] = Gℵ + Guℵ
end