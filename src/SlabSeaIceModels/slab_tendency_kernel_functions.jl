using ClimaSeaIce: latent_heat
using Oceananigans.Advection

# Thickness change due to accretion and melting, restricted by minimum allowable value
function ice_thickness_tendency(i, j, grid, clock,
                                velocities,
                                advection,
                                ice_thickness,
                                ice_consolidation_thickness,
                                top_temperature,
                                bottom_heat_bc,
                                top_external_heat_flux,
                                internal_heat_flux,
                                bottom_external_heat_flux,
                                phase_transitions,
                                h_forcing,
                                model_fields)

    Gh_advection = - div_Uc(i, j, 1, grid, advection, velocities, ice_thickness)

    @inbounds begin
        hᶜ = ice_consolidation_thickness[i, j, 1]
        hᵢ = ice_thickness[i, j, 1]
        Tuᵢ = top_temperature[i, j, 1]
    end

    # Consolidation criteria
    consolidated_ice = hᵢ >= hᶜ

    liquidus = phase_transitions.liquidus
    Tbᵢ = bottom_temperature(i, j, grid, bottom_heat_bc, liquidus)
    ℰb = latent_heat(phase_transitions, Tbᵢ)
    ℰu = latent_heat(phase_transitions, Tuᵢ)

    Qi = internal_heat_flux
    Qb = bottom_external_heat_flux
    Qu = top_external_heat_flux

    # Retrieve fluxes
    Quᵢ = getflux(Qu, i, j, grid, Tuᵢ, clock, model_fields)
    Qiᵢ = getflux(Qi, i, j, grid, Tuᵢ, clock, model_fields)
    Qbᵢ = getflux(Qb, i, j, grid, Tuᵢ, clock, model_fields)

    Qiᵢ = 0.0
    #Qbᵢ = 0.0

    #@show ℰb
    #@show ℰu

    # Compute forcing
    Fh = zero(grid) #h_forcing(i, j, grid, clock, model_fields)

    # If ice is consolidated, compute tendency for an ice slab; otherwise
    # just add ocean fluxes from frazil ice formation or melting
    slushy_Gh = - Qbᵢ / ℰb + Fh 

    # Upper (top) and bottom interface velocities
    wu = (Quᵢ - Qiᵢ) / ℰu # < 0 => melting
    wb = (Qiᵢ - Qbᵢ) / ℰb # < 0 => freezing

    slabby_Gh = wu + wb + Fh

    #@show consolidated_ice
    #@show Quᵢ
    #@show Qiᵢ
    #@show Qbᵢ
    #@show slabby_Gh
    #@show slushy_Gh
    #@show Fh

    return Gh_advection + ifelse(consolidated_ice, slabby_Gh, slushy_Gh)
end

#=
function tracer_tendency(i, j, grid, clock,
                         c_forcing,
                         c,
                         ice_thickness,
                         top_tracer_concentration,
                         bottom_tracer_concentration,
                         model_fields)

    wu = top_interface_velocity(i, j, grid, Tu, Qi, Qu, phase_transitions)
    wb = bottom_interface_velocity(i, j, grid, Tu, Qi, Qb, phase_transitions)

    cu = top_tracer_concentration[i, j]
    cb = bottom_tracer_concentration[i, j]

    return cu * wu - cb * wb + c_forcing(i, j, grid, clock, model_fields)
end
=#
