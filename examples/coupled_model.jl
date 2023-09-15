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

    grid = ocean.model.grid
    ice_ocean_thermal_flux = Field{Center, Center, Nothing}(grid)
    ice_ocean_salt_flux = Field{Center, Center, Nothing}(grid)

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
        # Ocean surface temperature
        T₀ = T[1, 1, Nz]
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

    # Initialization
    if coupled_model.ocean.model.clock.iteration == 0
        h⁻ = coupled_model.previous_ice_thickness
        hⁿ = coupled_model.ice.model.ice_thickness
        parent(h⁻) .= parent(hⁿ)
    end

    # TODO: put this in update_state!
    compute_ice_ocean_salinity_flux!(coupled_model)
    ice_ocean_latent_heat!(coupled_model)
    compute_solar_insolation!(coupled_model)
    compute_air_sea_flux!(coupled_model)

    time_step!(ice)
    time_step!(ocean)

    # TODO:
    # - Store fractional ice-free / ice-covered _time_ for more
    #   accurate flux computation?
    # - Or, input "excess heat flux" into ocean after the ice melts
    # - Currently, non-conservative for heat due bc we don't account for excess
        
    # TODO after ice time-step:
    #   - Adjust ocean temperature if the ice completely melts?
    
    return nothing
end

function compute_ice_ocean_salinity_flux!(coupled_model)
    # Compute salinity increment due to changes in ice thickness
    h⁻ = coupled_model.previous_ice_thickness
    hⁿ = coupled_model.ice.model.ice_thickness
    Sᵢ = coupled_model.ice.model.ice_salinity

    ocean = coupled_model.ocean
    Sₒ = ocean.model.tracers.S
    Qˢ = ocean.model.tracers.S.boundary_conditions.top.condition

    Nz = size(ocean.model.grid, 3)

    i = j = 1
    @inbounds begin
        # Thickness of surface grid cell
        Δz = Δzᶜᶜᶜ(i, j, Nz, ocean.model.grid)
        Δh = hⁿ[i, j, 1] - h⁻[i, j, 1]

        # Update surface salinity flux.
        # Note: the Δt below is the ocean time-step, eg.
        # ΔS = ⋯ - ∮ Qˢ dt ≈ ⋯ - Δtₒ * Qˢ 
        Qˢ[i, j, 1] = Δh / Δt * (Sᵢ[i, j, 1] - Sₒ[i, j, Nz])

        # Update previous ice thickness
        h⁻[i, j, 1] = hⁿ[i, j, 1]
    end

    return nothing
end

function ice_ocean_latent_heat!(coupled_model)
    ocean = coupled_model.ocean
    ice = coupled_model.ice
    ρₒ = coupled_model.ocean_density
    cₒ = coupled_model.ocean_heat_capacity
    Qₒ = ice.model.external_thermal_fluxes.bottom
    Tₒ = ocean.model.tracers.T
    Sₒ = ocean.model.tracers.S
    Δt = ocean.Δt
    hᵢ = ice.model.ice_thickness

    liquidus = ice.model.phase_transitions.liquidus
    grid = ocean.model.grid

    δQ = zero(grid)
    Nz = size(grid, 3)

    i = j = 1
    ice_covered = @inbounds hᵢ[i, j, 1] > 0.02 # 2 cm
    for k = Nz:-1:1
        @inbounds begin
            # Compute melting temperature
            Sᴺ = Sₒ[i, j, k]
            Tₘ = melting_temperature(liquidus, Sᴺ)

            # Compute total latent heat (per unit area) and latent heat flux
            Tᴺ = Tₒ[i, j, k]

            δE = ρₒ * cₒ * (Tᴺ - Tₘ) # < 0 in freezing conditions
                                     # > 0 in melting conditions
                                     
            freezing = Tᴺ < Tₘ
            δE = ifelse(ice_covered | freezing, δE, zero(grid))
            Tₒ[i, j, k] = ifelse(freezing, Tₘ, Tᴺ)

            # Tₒ[i, j, k] = ifelse((k == Nz) & ice_covered, Tₘ, Tᴺ)

            Δz = Δzᶜᶜᶜ(i, j, k, grid)
            δQ += δE * Δz / Δt # < 0 (we are warming the ocean)
        end
    end

    # Store ice-ocean flux (ignoring positive values computed when the, which 
    @inbounds Qₒ[i, j, 1] = δQ

    return nothing
end
