using ClimaSeaIce: latent_heat
using Oceananigans.Advection

# Thickness change due to accretion and melting, restricted by minimum allowable value
function thickness_tendency(i, j, grid, clock,
                            velocities,
                            advection,
                            thickness,
                            concentration,
                            consolidation_thickness,
                            top_temperature,
                            bottom_heat_bc,
                            top_external_heat_flux,
                            internal_heat_flux,
                            bottom_external_heat_flux,
                            phase_transitions,
                            h_forcing,
                            model_fields)

    Gh_advection = - div_Uc(i, j, 1, grid, advection, velocities, thickness)

    @inbounds begin
        hᶜ  = consolidation_thickness[i, j, 1]
        hᵢ  = thickness[i, j, 1]
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

    # Compute forcing
    Fh = zero(grid) #h_forcing(i, j, grid, clock, model_fields)

    # If ice is consolidated, compute tendency for an ice slab; otherwise
    # just add ocean fluxes from frazil ice formation or melting
    slushy_Gh = - Qbᵢ / ℰb + Fh 

    # Upper (top) and bottom interface velocities
    wu = (Quᵢ - Qiᵢ) / ℰu # < 0 => melting
    wb = (Qiᵢ - Qbᵢ) / ℰb # < 0 => freezing

    slabby_Gh = wu + wb + Fh

    return Gh_advection + ifelse(consolidated_ice, slabby_Gh, slushy_Gh)
end

# Concentration changes only due to advection?
function concentration_tendency(i, j, grid, clock,
                                velocities,
                                advection,
                                thickness,
                                concentration,
                                consolidation_thickness,
                                top_temperature,
                                bottom_heat_bc,
                                top_external_heat_flux,
                                internal_heat_flux,
                                bottom_external_heat_flux,
                                phase_transitions,
                                a_forcing,
                                model_fields)

    Gℵ_advection = - div_Uc(i, j, 1, grid, advection, velocities, concentration)

    return Gℵ_advection 
end

# Advection of sea-ice tracers
function tracer_tendency(i, j, grid, clock,
                         velocities,
                         advection,
                         thickness,
                         concentration,
                         tracer,
                         consolidation_thickness,
                         top_temperature,
                         bottom_heat_bc,
                         top_external_heat_flux,
                         internal_heat_flux,
                         bottom_external_heat_flux,
                         phase_transitions,
                         model_fields)

    Gc_advection = - div_Uc(i, j, 1, grid, advection, velocities, tracers)

    return Gc_advection
end
