using Oceananigans.Operators

struct IceOceanModel{FT, I, O, PI, PC}
    ice :: I
    previous_ice_thickness :: PI
    previous_ice_concentration :: PC
    ocean :: O
    ocean_density :: FT
    ocean_heat_capacity :: FT
    ocean_emissivity :: FT
    stefan_boltzmann_constant :: FT
    reference_temperature :: FT
end

function IceOceanModel(ice, ocean)
    previous_ice_thickness = deepcopy(ice.model.ice_thickness)
    previous_ice_concentration = deepcopy(ice.model.ice_concentration)

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
                         previous_ice_concentration,
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
    Nx, Ny, Nz = size(ocean.model.grid)

    for i = 1:Nx, j = 1:Ny
        @inbounds begin
            # Ocean surface temperature
            T₀ = T[i, j, Nz]
            h = ice.model.ice_thickness[i, j, 1]
            α = ice.model.ice_concentration[i, j, 1]
            I₀ = solar_insolation[i, j, 1]
        end

        # Radiation model
        ϵ = 1 # ocean emissivity
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

        @inbounds Qᵀ[i, j, 1] = (1 - α) * ΣQᵀ
    end

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
    αⁿ = coupled_model.ice.model.ice_concentration
    α⁻ = coupled_model.previous_ice_concentration
    Sᵢ = coupled_model.ice.model.ice_salinity

    ocean = coupled_model.ocean
    Sₒ = ocean.model.tracers.S
    Qˢ = ocean.model.tracers.S.boundary_conditions.top.condition

    Nx, Ny, Nz = size(ocean.model.grid)

    for i = 1:Nx, j = 1:Ny
        @inbounds begin
            # Thickness of surface grid cell
            Δz = Δzᶜᶜᶜ(i, j, Nz, ocean.model.grid)
            Δh = hⁿ[i, j, 1] - h⁻[i, j, 1]

            # Update surface salinity flux.
            # Note: the Δt below is the ocean time-step, eg.
            # ΔS = ⋯ - ∮ Qˢ dt ≈ ⋯ - Δtₒ * Qˢ 
            Qˢ[i, j, 1] = 0 #Δh / Δt * (Sᵢ[i, j, 1] - Sₒ[i, j, Nz])

            # Update previous ice thickness
            h⁻[i, j, 1] = hⁿ[i, j, 1]
        end
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
    αᵢ = ice.model.ice_concentration

    liquidus = ice.model.phase_transitions.liquidus
    grid = ocean.model.grid

    δQ = zero(grid)
    Nx, Ny, Nz = size(grid)

    for i = 1:Nx, j = 1:Ny
        icy_cell = @inbounds hᵢ[i, j, 1] > 0 # make ice bath approximation then

        for k = Nz:-1:1
            @inbounds begin
                # Various quantities
                Δz = Δzᶜᶜᶜ(i, j, k, grid)
                Tᴺ = Tₒ[i, j, k]
                Sᴺ = Sₒ[i, j, k]
            end

            # Melting / freezing temperature at the surface of the ocean
            Tₘ = melting_temperature(liquidus, Sᴺ)
                                     
            # Conditions for non-zero ice-ocean flux:
            #   - the ocean is below the freezing temperature, causing formation of ice.
            freezing = Tᴺ < Tₘ 

            #   - We are at the surface and the cell is covered by ice.
            icy_surface_cell = (k == Nz) & icy_cell

            # When there is a non-zero ice-ocean flux, we will instantaneously adjust the
            # temperature of the grid cells accordingly.
            adjust_temperature = freezing | icy_surface_cell

            # Compute change in ocean thermal energy.
            #
            #   - When Tᴺ < Tₘ, we heat the ocean back to melting temperature by extracting heat from the ice,
            #     assuming that the heat flux (which is carried by nascent ice crystals called frazil ice) floats
            #     instantaneously to the surface.
            #
            #   - When Tᴺ > Tₘ and we are in a surface cell covered by ice, we assume equilibrium
            #     and cool the ocean by injecting excess heat into the ice.
            # 
            δEₒ = adjust_temperature * ρₒ * cₒ * (Tₘ - Tᴺ)

            # Perform temperature adjustment
            @inline Tₒ[i, j, k] = ifelse(adjust_temperature, Tₘ, Tᴺ)

            # Compute the heat flux from ocean into ice.
            #
            # A positive value δQ > 0 implies that the ocean is cooled; ie heat
            # is fluxing upwards, into the ice. This occurs when applying the
            # ice bath equilibrium condition to cool down a warm ocean (δEₒ < 0).
            #
            # A negative value δQ < 0 implies that heat is fluxed from the ice into
            # the ocean, cooling the ice and heating the ocean (δEₒ > 0). This occurs when
            # frazil ice is formed within the ocean.
            
            δQ -= δEₒ * Δz / Δt
        end

        # Store ice-ocean flux
        @inbounds Qₒ[i, j, 1] = δQ
    end

    return nothing
end
