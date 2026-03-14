using Oceananigans.Architectures: architecture
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.Units: Time
using Oceananigans.Utils
using KernelAbstractions: @kernel, @index

thermodynamic_time_step!(model, ::Nothing, Δt) = nothing

function thermodynamic_time_step!(model, ::SlabThermodynamics, Δt)
    grid = model.grid
    arch = architecture(grid)

    if isnothing(model.snow_thermodynamics)
        launch!(arch, grid, :xy,
                _ice_thermodynamic_time_step!,
                model.ice_thickness,
                model.ice_concentration,
                grid, Δt,
                model.clock,
                model.ice_consolidation_thickness,
                model.ice_thermodynamics,
                model.external_heat_fluxes.top,
                model.external_heat_fluxes.bottom,
                fields(model))
    else
        launch!(arch, grid, :xy,
                _layered_thermodynamic_time_step!,
                model.ice_thickness,
                model.ice_concentration,
                grid, Δt,
                model.clock,
                model.ice_consolidation_thickness,
                model.ice_thermodynamics,
                model.external_heat_fluxes.top,
                model.external_heat_fluxes.bottom,
                model.snow_thickness,
                model.snow_thermodynamics,
                model.snow_precipitation,
                fields(model))
    end

    return nothing
end

#####
##### Ice-only thermodynamic kernel (no snow)
#####

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
@kernel function _ice_thermodynamic_time_step!(ice_thickness,
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

    ∂t_V = thermodynamic_tendency(i, j, 1, grid,
                                  ice_thermodynamics,
                                  ice_thickness,
                                  ice_concentration,
                                  ice_consolidation_thickness,
                                  top_external_heat_flux,
                                  bottom_external_heat_flux,
                                  clock, model_fields)

    hⁿ⁺¹, ℵⁿ⁺¹ = ice_volume_update(ice_thermodynamics, ∂t_V, hⁿ, ℵⁿ, hᶜ, Δt)

    @inbounds ice_concentration[i, j, 1] = ℵⁿ⁺¹
    @inbounds ice_thickness[i, j, 1]     = hⁿ⁺¹
end

#####
##### Layered snow + ice thermodynamic kernel
#####

# When snow is present, it sits on top of ice as an independent layer.
# The snow layer owns the surface temperature solve using the combined
# snow+ice conductive flux (IceSnowConductiveFlux). The interface temperature
# Tsi is computed analytically and written to the ice's top_surface_temperature.
#
# Top-down coupling: snow → Tsi → ice. The ice thermodynamics is unchanged.
@kernel function _layered_thermodynamic_time_step!(ice_thickness,
                                                    ice_concentration,
                                                    grid,
                                                    Δt,
                                                    clock,
                                                    ice_consolidation_thickness,
                                                    ice_thermodynamics,
                                                    top_external_heat_flux,
                                                    bottom_external_heat_flux,
                                                    snow_thickness,
                                                    snow_thermodynamics,
                                                    snow_precipitation,
                                                    model_fields)

    i, j = @index(Global, NTuple)

    @inbounds hⁿ  = ice_thickness[i, j, 1]
    @inbounds ℵⁿ  = ice_concentration[i, j, 1]
    @inbounds hᶜ  = ice_consolidation_thickness[i, j, 1]
    @inbounds hsⁿ = snow_thickness[i, j, 1]

    consolidated_ice = hⁿ ≥ hᶜ

    phase_ice = ice_thermodynamics.phase_transitions
    liquidus  = phase_ice.liquidus
    bottom_bc = ice_thermodynamics.heat_boundary_conditions.bottom
    @inbounds Si = model_fields.S[i, j, 1]

    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    Tm = melting_temperature(liquidus, Si)

    # --- Snow surface solve ---
    # Snow's internal_heat_flux is a FluxFunction wrapping ice_snow_conductive_flux,
    # which computes the combined conductive flux F = (Tb - Tu) / (hs/ks + hi/ki).
    # The MeltingConstrainedFluxBalance root-find balances Qatm(Tu) = Fc(Tu).
    snow_top_bc = snow_thermodynamics.heat_boundary_conditions.top
    Tu = snow_thermodynamics.top_surface_temperature
    Qi = snow_thermodynamics.internal_heat_flux
    Qu = top_external_heat_flux

    # Effective melting temperature: snow (0 C) when snow present, ice otherwise
    Tm = ifelse(hsⁿ > 0, zero(Tm), Tm)

    if !isa(snow_top_bc, PrescribedTemperature)
        if consolidated_ice
            Tu⁻ = @inbounds Tu[i, j, 1]
            Tuⁿ = top_surface_temperature(i, j, grid, snow_top_bc, Tu⁻, Qi, Qu, clock, model_fields)
            Tuⁿ = min(Tuⁿ, Tm)
        else
            Tuⁿ = Tb
        end
        @inbounds Tu[i, j, 1] = Tuⁿ
    end

    @inbounds Tus = Tu[i, j, 1]

    # Tsi = Tb + (Tus - Tb) * Ri / (Rs + Ri)
    # When hs = 0: Tsi = Tus (snow layer has zero resistance)
    ks  = snow_thermodynamics.internal_heat_flux.parameters.snow_conductivity
    Tsi = interface_temperature(i, j, grid, ice_thermodynamics, ks, Tus, model_fields)
    @inbounds ice_thermodynamics.top_surface_temperature[i, j, 1] = Tsi

    Qis = ifelse(consolidated_ice, getflux(Qi, i, j, grid, Tus, clock, model_fields), zero(grid))
    Qui = getflux(Qu, i, j, grid, Tus, clock, model_fields)

    # δQ < 0 means more internal flux than external → energy available for melting
    δQ = Qui - Qis
    melt_energy = max(zero(δQ), -δQ) # positive when melting

    ρs = snow_thermodynamics.phase_transitions.density
    Ls = snow_thermodynamics.phase_transitions.reference_latent_heat

    snow_energy_capacity = ρs * Ls * hsⁿ / Δt # W/m²
    Qs  = min(melt_energy, snow_energy_capacity)
    Gs⁻ = Qs / (ρs * Ls)

    # The effective top flux for the ice is Qui + snow_absorbed.
    # Snow absorbing melt energy acts as extra cooling from the ice's perspective:
    # thermodynamic_tendency computes Qi_ice = Fc (by interface temperature construction),
    # so wu = (Qui + snow_absorbed - Fc) / ℰu = -excess / ℰu
    # and wb = (Fc - Qb) / ℰb.
    Qui = Qui + Qs

    ∂t_V = thermodynamic_tendency(i, j, 1, grid,
                                  ice_thermodynamics,
                                  ice_thickness,
                                  ice_concentration,
                                  ice_consolidation_thickness,
                                  Qui,
                                  bottom_external_heat_flux,
                                  clock, model_fields)

    hⁿ⁺¹, ℵⁿ⁺¹ = ice_volume_update(ice_thermodynamics, ∂t_V, hⁿ, ℵⁿ, hᶜ, Δt)

    # Conserve snow volume when concentration changes: new ice has no snow,
    # so hs adjusts to keep hs * ℵ constant (analogous to how ice tracks h * ℵ).
    hsⁿ = ifelse(ℵⁿ⁺¹ > 0, hsⁿ * ℵⁿ / ℵⁿ⁺¹, zero(hsⁿ))

    Gs⁺ = snow_accumulation(i, j, snow_precipitation, snow_thermodynamics, ℵⁿ⁺¹, clock)
    hs⁺  = hsⁿ + Δt * (Gs⁺ - Gs⁻)
    hs⁺  = max(zero(hs⁺), hs⁺)

    # Snow-ice formation (flooding when freeboard is negative)
    hⁿ⁺¹, hs⁺ = snow_ice_formation(hⁿ⁺¹, hs⁺, ice_thermodynamics, snow_thermodynamics)

    # Reset snow when no ice
    hs⁺ = ifelse(ℵⁿ⁺¹ ≤ 0, zero(hs⁺), hs⁺)

    @inbounds ice_concentration[i, j, 1] = ℵⁿ⁺¹
    @inbounds ice_thickness[i, j, 1]     = hⁿ⁺¹
    @inbounds snow_thickness[i, j, 1]    = hs⁺
end

#####
##### Shared helper functions
#####

@inline function ice_volume_update(ice_thermodynamics, ∂t_V, hⁿ, ℵⁿ, hᶜ, Δt)
    Vⁿ⁺¹ = hⁿ * ℵⁿ + Δt * ∂t_V
    Vⁿ⁺¹ = max(zero(Vⁿ⁺¹), Vⁿ⁺¹)

    ∂t_V = (Vⁿ⁺¹ - hⁿ * ℵⁿ) / Δt
    ℵ⁺   = concentration_thermodynamic_step(ice_thermodynamics.concentration_evolution, ∂t_V, ℵⁿ, hⁿ, hᶜ, Δt)
    h⁺   = Vⁿ⁺¹ / ℵ⁺

    # Treat pathological cases
    h⁺ = ifelse(ℵ⁺ ≤ 0, zero(h⁺), h⁺)
    ℵ⁺ = ifelse(∂t_V == 0, ℵⁿ, ℵ⁺)     # No volume change
    h⁺ = ifelse(∂t_V == 0, hⁿ, h⁺)     # No volume change
    ℵ⁺ = ifelse(h⁺ == 0, zero(ℵ⁺), ℵ⁺) # reset the concentration if there is no sea-ice
    h⁺ = ifelse(ℵ⁺ == 0, zero(h⁺), h⁺) # reset the thickness if there is no sea-ice

    # Ridging caused by the thermodynamic step
    ℵⁿ⁺¹ = ifelse(ℵ⁺ > 1, one(ℵ⁺), ℵ⁺)
    hⁿ⁺¹ = ifelse(ℵ⁺ > 1,  h⁺ * ℵ⁺, h⁺)

    return (hⁿ⁺¹, ℵⁿ⁺¹)
end

@inline get_precipitation(i, j, snow_precip, clock) = @inbounds snow_precip[i, j, 1]
@inline get_precipitation(i, j, snow_precip::FieldTimeSeries, clock) = @inbounds snow_precip[i, j, 1, Time(clock.time)]

@inline function snow_accumulation(i, j, snow_precip, snow_thermo, ℵ, clock)
    Ps = get_precipitation(i, j, snow_precip, clock) # kg/m^2/s
    ρs = snow_thermo.phase_transitions.density
    return ifelse(ℵ > 0, Ps / ρs, zero(ρs))
end

@inline function snow_ice_formation(hi, hs, ice_thermo, snow_thermo)
    ρi = ice_thermo.phase_transitions.density
    ρs = snow_thermo.phase_transitions.density
    ρw = ice_thermo.phase_transitions.liquid_density

    # Freeboard: positive when ice floats above waterline
    hf = hi * (1 - ρi / ρw) - hs * ρs / ρw

    # Energy-conserving snow-ice formation.
    # Since E = -ℵ(ℰ h + ρs Ls hs) with ℰ = ρw L (see `latent_heat`),
    # energy conservation requires ℰ δhi = ρs Ls δhs, i.e. δhs = ρw δhi / ρs
    # (with L = Ls). Combined with freeboard → 0:
    #   δhi = -hf / (2 - ρi/ρw)
    flooding = hf < 0
    denom = 2 - ρi / ρw
    δhi = ifelse(flooding, -hf / denom, zero(hf))
    δhs = ifelse(flooding, ρw * δhi / ρs, zero(hf))

    # Clip: can't remove more snow than available
    hs⁺ = max(zero(hs), hs - δhs)
    δhs = hs - hs⁺

    # Recompute δhi from actual snow removed (energy conservation)
    δhi = δhs * ρs / ρw
    hi⁺ = hi + δhi

    return (hi⁺, hs⁺)
end

# We parameterize the evolution of ice thickness and concentration
# (i.e. lateral vs vertical growth) following Hibler (1979)
@inline function concentration_thermodynamic_step(::ProportionalEvolution, ∂t_V, ℵⁿ, hⁿ, hᶜ, Δt)
    freezing = (∂t_V ≥ 0) # Freezing
    melting  = (∂t_V < 0) # Melting

    ∂t_ℵᶠ = (1 - ℵⁿ) /  hᶜ * ∂t_V * freezing
    ∂t_ℵᵐ =      ℵⁿ  / 2hⁿ * ∂t_V * melting

    # Update concentration accordingly
    ℵ⁺ = ℵⁿ + Δt * (∂t_ℵᶠ + ∂t_ℵᵐ)
    ℵ⁺ = max(zero(ℵ⁺), ℵ⁺)

    return ℵ⁺
end
