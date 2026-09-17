using KernelAbstractions: @kernel, @index
using Oceananigans.Architectures: architecture
using Oceananigans.OutputReaders: FieldTimeSeries, GPUAdaptedFieldTimeSeries
using Oceananigans.Units: Time

# No ice, no thermodynamics!
thermodynamic_time_step!(model, ::Nothing, snow_thermodynamics, Δt) = nothing

# Only slab ice
function thermodynamic_time_step!(model, ice_thermodynamics::SlabThermodynamics, ::Nothing, Δt)
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy,
            _ice_thermodynamic_time_step!,
            model.ice_thickness,
            model.ice_concentration,
            grid, Δt,
            model.clock,
            model.ice_consolidation_thickness,
            ice_thermodynamics,
            model.phase_transitions,
            model.sea_ice_density,
            model.external_heat_fluxes.top,
            model.external_heat_fluxes.bottom,
            model.mass_fluxes,
            fields(model))

    return nothing
end

# Slab ice and slab snow
function thermodynamic_time_step!(model,
                                  ice_thermodynamics::SlabThermodynamics,
                                  snow_thermodynamics::SlabThermodynamics,
                                  Δt)
    grid = model.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy,
            _layered_thermodynamic_time_step!,
            model.ice_thickness,
            model.ice_concentration,
            grid, Δt,
            model.clock,
            model.ice_consolidation_thickness,
            ice_thermodynamics,
            snow_thermodynamics,
            model.phase_transitions,
            model.sea_ice_density,
            model.snow_density,
            model.external_heat_fluxes.top,
            model.external_heat_fluxes.bottom,
            model.snow_thickness,
            model.snowfall,
            model.mass_fluxes,
            fields(model))

    return nothing
end

#####
##### Ice-only thermodynamic kernel (no snow)
#####

# The thermodynamic step is computed in a single kernel following:
#
#     ∂t_V = ∂t_h * ℵ + h * ∂t_ℵ
#
# Therefore:
#
#     ∂t_V = Gᴸ + Gⱽ = (h⁺ * ℵ⁺ - hⁿ * ℵⁿ) / Δt
#
# The two will be adjusted conservatively after the thermodynamic step to ensure that ℵ ≤ 1.
@kernel function _ice_thermodynamic_time_step!(ice_thickness,
                                               ice_concentration,
                                               grid,
                                               Δt,
                                               clock,
                                               ice_consolidation_thickness,
                                               ice_thermodynamics,
                                               phase_transitions,
                                               sea_ice_density,
                                               top_external_heat_flux,
                                               bottom_external_heat_flux,
                                               mass_fluxes,
                                               model_fields)

    i, j = @index(Global, NTuple)

    @inbounds hⁿ = ice_thickness[i, j, 1]
    @inbounds ℵⁿ = ice_concentration[i, j, 1]
    @inbounds hᶜ = ice_consolidation_thickness[i, j, 1]

    ∂t_𝓋 = thermodynamic_tendency(i, j, 1, grid,
                                  ice_thermodynamics,
                                  phase_transitions,
                                  sea_ice_density,
                                  ice_thickness,
                                  ice_concentration,
                                  ice_consolidation_thickness,
                                  top_external_heat_flux,
                                  bottom_external_heat_flux,
                                  clock, model_fields)

    hⁿ⁺¹, ℵⁿ⁺¹ = ice_volume_update(ice_thermodynamics, ∂t_𝓋, hⁿ, ℵⁿ, hᶜ, Δt)

    @inbounds begin
        ρi = sea_ice_density[i, j, 1]
        
        ice_concentration[i, j, 1] = ℵⁿ⁺¹
        ice_thickness[i, j, 1]     = hⁿ⁺¹

        mass_fluxes.thermodynamics.ice[i, j, 1]  = ρi * (hⁿ⁺¹ * ℵⁿ⁺¹ - hⁿ * ℵⁿ) / Δt
        mass_fluxes.thermodynamics.snow[i, j, 1] = 0
        mass_fluxes.intercepted_snowfall[i, j, 1] = 0
    end
end

#####
##### Layered snow + ice thermodynamic kernel
#####

# When snow is present, it sits on top of ice as an independent layer. The
# column-top surface solve happens at the atmosphere-snow surface using the
# combined snow+ice conductive flux (`IceSnowConductiveFlux`).
#
# The snow-ice interface temperature Tsi is computed analytically from the
# resistance ratio and handed to `ice_melt_freeze_tendency` as an argument,
# so the ice's surface solve is not invoked.
@kernel function _layered_thermodynamic_time_step!(ice_thickness,
                                                   ice_concentration,
                                                   grid,
                                                   Δt,
                                                   clock,
                                                   ice_consolidation_thickness,
                                                   ice_thermodynamics,
                                                   snow_thermodynamics,
                                                   phase_transitions,
                                                   sea_ice_density,
                                                   snow_density,
                                                   top_external_heat_flux,
                                                   bottom_external_heat_flux,
                                                   snow_thickness,
                                                   snowfall,
                                                   mass_fluxes,
                                                   model_fields)

    i, j = @index(Global, NTuple)

    @inbounds hiⁿ = ice_thickness[i, j, 1]
    @inbounds ℵⁿ  = ice_concentration[i, j, 1]
    @inbounds hᶜ  = ice_consolidation_thickness[i, j, 1]
    @inbounds hsⁿ = snow_thickness[i, j, 1]

    # Per-cell volumes before the step, captured before `hsⁿ` is rebased below
    Viⁿ = hiⁿ * ℵⁿ
    Vsⁿ = hsⁿ * ℵⁿ

    consolidated_ice = hiⁿ ≥ hᶜ

    liquidus = phase_transitions.liquidus
    bottom_bc = ice_thermodynamics.heat_boundary_conditions.bottom
    @inbounds Si = model_fields.S[i, j, 1]

    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    Tm = melting_temperature(liquidus, Si)

    # Snow surface solve using the combined snow+ice conductive flux with a
    # resistors in series formulation: F = (Tb - Tu) / (hs/ks + hi/ki).
    ks = snow_thermodynamics.internal_heat_flux.conductivity
    ki = ice_thermodynamics.internal_heat_flux.conductivity
    combined_flux = IceSnowConductiveFlux(ks, ki)

    # Column internal heat flux
    Qic = internal_flux_function(combined_flux, liquidus, bottom_bc)

    # Ice-only flux wrapper for the ice-interior evaluation at Tsi.
    Qii = internal_flux_function(ice_thermodynamics.internal_heat_flux, liquidus, bottom_bc)

    snow_top_bc = snow_thermodynamics.heat_boundary_conditions.top
    Tu = snow_thermodynamics.top_surface_temperature
    Qu = top_external_heat_flux

    # Effective melting temperature: snow (0 C) when snow present, ice otherwise
    Tm = ifelse(hsⁿ > 0, zero(Tm), Tm)

    if !isa(snow_top_bc, PrescribedTemperature)
        if consolidated_ice
            @inbounds Tu⁻ = Tu[i, j, 1]
            Tuⁿ = top_surface_temperature(i, j, grid, snow_top_bc, Tu⁻, Qic, Qu, clock, model_fields)
            Tuⁿ = min(Tuⁿ, Tm)
        else
            Tuⁿ = Tb
        end
        @inbounds Tu[i, j, 1] = Tuⁿ
    end

    @inbounds Tus = Tu[i, j, 1]

    # Tsi = Tb + (Tus - Tb) * Ri / (Rs + Ri)
    # When hs = 0: Tsi = Tus (snow layer has zero resistance)
    Tsi = interface_temperature(i, j, grid, combined_flux, bottom_bc, liquidus, Tus, model_fields)

    # Store Tsi as the ice top surface temperature
    @inbounds ice_thermodynamics.top_surface_temperature[i, j, 1] = Tsi

    # Snow-surface energy balance (Qis per-ice column flux, Qui per-cell from
    # the coupler). Converting Qui to per-ice...
    Qis = ifelse(consolidated_ice, getflux(Qic, i, j, grid, Tus, clock, model_fields), zero(grid))
    Qui = getflux(Qu, i, j, grid, Tus, clock, model_fields)

    Qui_per_ice = ifelse(ℵⁿ > 0, Qui / ℵⁿ, zero(Qui))

    δQ = Qui_per_ice - Qis                       # δQ < 0 ⇒ energy available for melt
    melt_energy = max(zero(δQ), -δQ)             # per-ice

    @inbounds ρs = snow_density[i, j, 1]
    ℒs = phase_transitions.reference_latent_heat

    snow_energy_capacity = ρs * ℒs * hsⁿ / Δt    # per-ice
    Qs  = min(melt_energy, snow_energy_capacity) # per-ice
    Gs⁻ = Qs / (ρs * ℒs)                         # per-ice, drives Δhs

    # Closed-form self-consistent solve for ℵⁿ⁺¹. The ice-top effective flux
    # Quiᵉᶠᶠ = Qui + Qs·ℵⁿ⁺¹ couples Quiᵉᶠᶠ and ℵⁿ⁺¹ linearly; the slab's
    # concentration rule is also linear in ∂t_V, so the fixed point
    #   ∂t_V = α + β·ℵⁿ⁺¹,    α = (Qui − Qbi) / (ρᵢℒ),  β = Qs/(ρᵢℒ)
    #   ℵⁿ⁺¹ = ℵⁿ + K·∂t_V,   K = Δt · C,
    #   C = ℵⁿ/(2hⁿ) (melt)  or  (1−ℵⁿ)/hᶜ (freeze)
    # has the explicit solution  ℵⁿ⁺¹ = (ℵⁿ + K·α) / (1 − K·β).
    # Solve both branches and pick the one with sign(∂t_𝓋(ℵⁿ⁺¹)) consistent
    # with its branch assumption.
    @inbounds ρi = sea_ice_density[i, j, 1]
    ρiℒ = ρi * ℒs
    Qbi = getflux(bottom_external_heat_flux, i, j, grid, Tus, clock, model_fields)

    α = (Qui - Qbi) / ρiℒ  # ice-content-per-area rate (∂t_𝓋, m/s)
    β = Qs / ρiℒ           # coefficient on ℵⁿ⁺¹

    Cᵐ = ifelse(hiⁿ > zero(hiⁿ), ℵⁿ / (2 * hiⁿ), zero(hiⁿ))
    Cᶠ = ifelse(hᶜ  > zero(hᶜ),  (1 - ℵⁿ) / hᶜ,  zero(hᶜ))
    Kᵐ = Δt * Cᵐ
    Kᶠ = Δt * Cᶠ

    ε  = eps(typeof(β))
    Dᵐ = 1 - Kᵐ * β
    Dᶠ = 1 - Kᶠ * β
    ℵᵐ = ifelse(abs(Dᵐ) > ε, (ℵⁿ + Kᵐ * α) / Dᵐ, ℵⁿ + Kᵐ * α)
    ℵᶠ = ifelse(abs(Dᶠ) > ε, (ℵⁿ + Kᶠ * α) / Dᶠ, ℵⁿ + Kᶠ * α)

    # Branch selection: melt if melt-branch solution yields ∂t_𝓋 < 0, else freeze.
    ∂t_𝓋ᵐ   = α + β * ℵᵐ
    melting = ∂t_𝓋ᵐ < zero(∂t_𝓋ᵐ)
    ℵtmp    = ifelse(melting, ℵᵐ, ℵᶠ)

    # Final state via `ice_volume_update` (handles 𝓋<0 clipping, ridging, etc.).
    # Pass the cached Qbi scalar (not the closure) so the bottom-flux closure is evaluated exactly once per step
    Quiᵉᶠᶠ = Qui + Qs * ℵtmp
    ∂t_𝓋 = ice_melt_freeze_tendency(i, j, 1, grid,
                                    ice_thermodynamics,
                                    phase_transitions,
                                    sea_ice_density,
                                    Qii,
                                    Tsi,
                                    ice_thickness, ice_consolidation_thickness,
                                    Quiᵉᶠᶠ, Qbi,
                                    clock, model_fields)

    hiⁿ⁺¹, ℵⁿ⁺¹ = ice_volume_update(ice_thermodynamics, ∂t_𝓋, hiⁿ, ℵⁿ, hᶜ, Δt)

    # Conserve snow volume when concentration changes: new ice has no snow,
    # so hs adjusts to keep hs * ℵ constant (analogous to how ice tracks hi * ℵ).
    hsⁿ = ifelse(ℵⁿ⁺¹ > 0, hsⁿ * ℵⁿ / ℵⁿ⁺¹, zero(hsⁿ))

    Gs⁺ = snow_accumulation(i, j, snowfall, ρs, ℵⁿ⁺¹, clock)
    hs⁺ = hsⁿ + Δt * (Gs⁺ - Gs⁻)
    hs⁺ = max(zero(hs⁺), hs⁺)

    # Snow-ice formation (flooding when freeboard is negative).
    hiⁿ⁺¹, hs⁺ = snow_ice_formation(hiⁿ⁺¹, hs⁺, ρi, ρs, phase_transitions.liquid_density)

    # Reset snow when no ice
    hs⁺ = ifelse(ℵⁿ⁺¹ ≤ 0, zero(hs⁺), hs⁺)

    @inbounds ice_concentration[i, j, 1] = ℵⁿ⁺¹
    @inbounds ice_thickness[i, j, 1]     = hiⁿ⁺¹
    @inbounds snow_thickness[i, j, 1]    = hs⁺

    # Snowfall mass deposited on top of the ice
    Psᵃᵇˢ = ρs * Gs⁺ * ℵⁿ⁺¹

    @inbounds begin
        mass_fluxes.thermodynamics.ice[i, j, 1]  = ρi * (hiⁿ⁺¹ * ℵⁿ⁺¹ - Viⁿ) / Δt
        mass_fluxes.thermodynamics.snow[i, j, 1] = ρs * (hs⁺ * ℵⁿ⁺¹ - Vsⁿ) / Δt - Psᵃᵇˢ
        mass_fluxes.intercepted_snowfall[i, j, 1] = Psᵃᵇˢ
    end
end

#####
##### Shared helper functions
#####

@inline function ice_volume_update(ice_thermodynamics, ∂t_𝓋, hⁿ, ℵⁿ, hᶜ, Δt)
    𝓋ⁿ⁺¹ = hⁿ * ℵⁿ + Δt * ∂t_𝓋
    𝓋ⁿ⁺¹ = max(zero(𝓋ⁿ⁺¹), 𝓋ⁿ⁺¹)

    ∂t_𝓋 = (𝓋ⁿ⁺¹ - hⁿ * ℵⁿ) / Δt
    ℵ⁺   = concentration_thermodynamic_step(ice_thermodynamics.concentration_evolution, ∂t_𝓋, ℵⁿ, hⁿ, hᶜ, Δt)
    h⁺   = 𝓋ⁿ⁺¹ / ℵ⁺

    # Treat pathological cases
    h⁺ = ifelse(ℵ⁺ ≤ 0, zero(h⁺), h⁺)
    ℵ⁺ = ifelse(∂t_𝓋 == 0, ℵⁿ, ℵ⁺)     # No volume change
    h⁺ = ifelse(∂t_𝓋 == 0, hⁿ, h⁺)     # No volume change
    ℵ⁺ = ifelse(h⁺ == 0, zero(ℵ⁺), ℵ⁺) # reset the concentration if there is no sea-ice
    h⁺ = ifelse(ℵ⁺ == 0, zero(h⁺), h⁺) # reset the thickness if there is no sea-ice

    # Ridging caused by the thermodynamic step
    ℵⁿ⁺¹ = ifelse(ℵ⁺ > 1, one(ℵ⁺), ℵ⁺)
    hⁿ⁺¹ = ifelse(ℵ⁺ > 1,  h⁺ * ℵ⁺, h⁺)

    return (hⁿ⁺¹, ℵⁿ⁺¹)
end

const FTS = Union{FieldTimeSeries, GPUAdaptedFieldTimeSeries}

@inline get_precipitation(i, j, Ps,      clock) = @inbounds Ps[i, j, 1]
@inline get_precipitation(i, j, Ps::FTS, clock) = @inbounds Ps[i, j, 1, Time(clock.time)]

@inline function snow_accumulation(i, j, snowfall, ρs, ℵ, clock)
    Ps = get_precipitation(i, j, snowfall, clock) # kg m⁻² s⁻¹
    return ifelse(ℵ > 0, Ps / ρs, zero(ρs))
end

@inline function snow_ice_formation(hi, hs, ρi, ρs, ρw)
    # Freeboard: positive when ice floats above waterline
    hf = hi * (1 - ρi / ρw) - hs * ρs / ρw

    # Flooding occurs when freeboard is negative.
    # Mass conservation (δhi ρi = δhs ρs) also conserves energy
    # since the per-mass latent heat is shared between snow and ice.
    flooding = hf < 0
    δhs = ifelse(flooding, -hf * ρi / ρs, zero(hf))

    # Clip: can't remove more snow than available
    hs⁺ = max(zero(hs), hs - δhs)
    δhs = hs - hs⁺
    # Recompute δhi from actual snow removed (energy conservation)
    δhi = δhs * ρs / ρi
    hi⁺ = hi + δhi

    return (hi⁺, hs⁺)
end

# We parameterize the evolution of ice thickness and concentration
# (i.e. lateral vs vertical growth) following Hibler (1979)
@inline function concentration_thermodynamic_step(::ProportionalEvolution, ∂t_𝓋, ℵⁿ, hⁿ, hᶜ, Δt)
    freezing = (∂t_𝓋 ≥ 0) # Freezing
    melting  = (∂t_𝓋 < 0) # Melting

    ∂t_ℵᶠ = (1 - ℵⁿ) /  hᶜ * ∂t_𝓋 * freezing
    ∂t_ℵᵐ =      ℵⁿ  / 2hⁿ * ∂t_𝓋 * melting

    # Update concentration accordingly
    ℵ⁺ = ℵⁿ + Δt * (∂t_ℵᶠ + ∂t_ℵᵐ)
    ℵ⁺ = max(zero(ℵ⁺), ℵ⁺)

    return ℵ⁺
end
