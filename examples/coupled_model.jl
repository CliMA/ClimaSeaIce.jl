using Oceananigans.Operators

struct IceOceanModel{FT, I, O, PI}
    ice :: I
    previous_ice_thickness:: PI
    ocean :: O
    ocean_density :: FT
    ocean_heat_capacity :: FT
    ocean_emissivity :: FT
    stefan_boltzmann_constant :: FT
    reference_temperature :: FT
end

function IceOceanModel(ice, ocean)
    previous_ice_thickness = deepcopy(ice.model.ice_thickness)

    ocean_density = 1024
    ocean_heat_capacity = 3991
    ocean_emissivity = 1
    reference_temperature = 273.15
    stefan_boltzmann_constant = 5.67e-8

    # How would we ensure consistency?
    try
        if ice.model.external_thermal_fluxes.top isa RadiativeEmission
            radiation = ice.model.external_thermal_fluxes.top
        else
            radiation = filter(flux isa RadiativeEmission, ice.model.external_thermal_fluxes.top) |> first
        end

        stefan_boltzmann_constant = radiation.stefan_boltzmann_constant
        reference_temperature = radiation.reference_temperature
    catch
    end

    FT = eltype(ocean.model.grid)

    return IceOceanModel(ice,
                         previous_ice_thickness,
                         ocean,
                         convert(FT, ocean_density),
                         convert(FT, ocean_heat_capacity),
                         convert(FT, ocean_emissivity),
                         convert(FT, stefan_boltzmann_constant),
                         convert(FT, reference_temperature))
end

time(coupled_model::IceOceanModel) = time(coupled_model.ocean)

function compute_air_sea_flux!(coupled_model)
    ocean = coupled_model.ocean
    ice = coupled_model.ice
    radiation = ice.model.external_thermal_fluxes.top[2]

    T = ocean.model.tracers.T
    Nz = size(ocean.model.grid, 3)

    @inbounds begin
        T₀ = T[1, 1, 1]
        h = ice.model.ice_thickness[1, 1, 1]
        I₀ = solar_insolation[1, 1, 1]
    end

    # Radiation model
    ϵ = 1.0 # ocean emissivity
    σ = coupled_model.stefan_boltzmann_constant
    ρₒ = coupled_model.ocean_density
    cₒ = coupled_model.ocean_heat_capacity
    Tᵣ = coupled_model.reference_temperature
    ΣQᵀ = ϵ * σ * (T₀ + Tᵣ)^4 / (ρₒ * cₒ)

    # Also add solar insolation
    ΣQᵀ += I₀ / (ρₒ * cₒ)

    # Set the surface flux only if ice-free
    Qᵀ = T.boundary_conditions.top.condition
    grid = ocean.model.grid
    @inbounds Qᵀ[1, 1, 1] = ifelse(h > 0, zero(grid), ΣQᵀ)

    return nothing
end

function time_step!(coupled_model::IceOceanModel, Δt)
    ocean = coupled_model.ocean
    ice = coupled_model.ice
    liquidus = ice.model.phase_transitions.liquidus
    grid = ocean.model.grid
    ice.Δt = Δt
    ocean.Δt = Δt

    compute_solar_insolation!(coupled_model)
    compute_air_sea_flux!(coupled_model)

    time_step!(ice)

    # TODO:
    # - Store fractional ice-free / ice-covered _time_ for more
    #   accurate flux computation?
    # - Or, input "excess heat flux" into ocean after the ice melts
    # - Currently, non-conservative for heat due bc we don't account for excess

    # Compute salinity increment due to changes in ice thickness
    h⁻ = coupled_model.previous_ice_thickness
    hⁿ = coupled_model.ice.model.ice_thickness
    Sᵢ = coupled_model.ice.model.ice_salinity
    Sₒ = ocean.model.tracers.S

    @inbounds begin
        Δz = Δzᶜᶜᶜ(1, 1, 1, ocean.model.grid)
        Δh = hⁿ[1, 1, 1] - h⁻[1, 1, 1]
        ΔS = Δh / Δz * (Sₒ[1, 1, 1] - Sᵢ[1, 1, 1])

        # Update ocean salinity
        Sₒ[1, 1, 1] += ΔS

        # Update previous ice thickness
        h⁻[1, 1, 1] = hⁿ[1, 1, 1]
    end

    time_step!(ocean)

    # Compute ice-ocean latent heat flux
    ρₒ = coupled_model.ocean_density
    cₒ = coupled_model.ocean_heat_capacity
    hₒ = ocean.model.grid.Lz # mixed layer depth
    Qₒ = ice.model.external_thermal_fluxes.bottom
    Tₒ = ocean.model.tracers.T

    @inbounds begin
        T₁ = Tₒ[1, 1, 1]
        S₁ = Sₒ[1, 1, 1]

        # Compute total latent heat (per unit area) and latent heat flux
        Tₘ = melting_temperature(liquidus, S₁)
        δE = ρₒ * cₒ * (T₁ - Tₘ) # > 0
        δQ = δE * hₒ / Δt        # > 0 (we are warming the ocean)

        # Clip 
        Qₒ[1, 1, 1] = min(zero(grid), δQ)
        Tₒ[1, 1, 1] = max(Tₘ, Tₒ[1, 1, 1])
    end

    # TODO after ice time-step:
    #   - Adjust ocean temperature if the ice completely melts?
    
    return nothing
end

