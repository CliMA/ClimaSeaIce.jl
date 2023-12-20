using ClimaSeaIce: latent_heat
using Oceananigans.Advection
using Oceananigans.Operators
using Oceananigans.Advection: _advective_tracer_flux_x, _advective_tracer_flux_y

@inline horizontal_div_Uc(i, j, grid, ::Nothing, U, c) = zero(grid)
@inline horizontal_div_Uc(i, j, grid, advection, U, c) = 
    1 / Vᶜᶜᶜ(i, j, 1, grid) * (δxᶜᵃᵃ(i, j, 1, grid, _advective_tracer_flux_x, advection, U.u, c) +
                               δyᵃᶜᵃ(i, j, 1, grid, _advective_tracer_flux_y, advection, U.v, c))
                               
@kernel function _compute_tracer_tendencies!(tendencies,
                                             grid,
                                             clock,
                                             thickness,
                                             concentration,
                                             velocities,
                                             advection,
                                             top_temperature,
                                             top_heat_bc,
                                             bottom_heat_bc,
                                             top_external_heat_flux,
                                             internal_heat_flux,
                                             bottom_external_heat_flux,
                                             consolidation_thickness,
                                             phase_transitions,
                                             forcing,
                                             model_fields)

    i, j = @index(Global, NTuple)

    ℵ  = concentration
    h  = thickness
    hᶜ = consolidation_thickness
    Qi = internal_heat_flux
    Qu = top_external_heat_flux
    Qb = bottom_external_heat_flux
    Tu = top_temperature
    liquidus = phase_transitions.liquidus
    
    Gh = tendencies.h
    Gℵ = tendencies.ℵ

    # Determine top surface temperature
    if !isa(top_heat_bc, PrescribedTemperature) # update surface temperature?

        consolidated_ice = @inbounds h[i, j, 1] >= hᶜ[i, j, 1]

        if consolidated_ice # slab is consolidated and has an independent surface temperature
            Tu⁻ = @inbounds Tu[i, j, 1]
            Tuⁿ = top_surface_temperature(i, j, grid, top_heat_bc, Tu⁻, Qi, Qu, clock, model_fields)
        else # slab is unconsolidated and does not have an independent surface temperature
            Tuⁿ = bottom_temperature(i, j, grid, bottom_heat_bc, liquidus)
        end

        @inbounds Tu[i, j, 1] = Tuⁿ
    end

    tracer_args = (i, j, grid, clock, velocities,
                   advection,
                   thickness,
                   concentration,
                   consolidation_thickness,
                   top_temperature,
                   bottom_heat_bc,
                   top_external_heat_flux,
                   internal_heat_flux,
                   bottom_external_heat_flux,
                   phase_transitions)

    @inbounds Gh[i, j, 1] = thickness_tendency(tracer_args..., nothing, model_fields)
    @inbounds Gℵ[i, j, 1] = concentration_tendency(tracer_args..., nothing, model_fields)
end

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

    Gh_advection = - div_Uℵh(i, j, grid, advection, velocities, concentration, thickness)

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

    Gℵ_advection = - horizontal_div_Uc(i, j, grid, advection, velocities, concentration)

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

    return  - div_Uh(i, j, grid, advection, velocities, concentration, tracer)
end
