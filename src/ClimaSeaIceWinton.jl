module ClimaSeaIceWinton

using ClimaSeaIce

# Do we want to create an abstract ThermodynamicSeaIceModel struct?
# For now I'm making this a subtype of ThermodynamicSeaIceModel
#mutable struct WintonThermodynamicSeaIceModel{} <: ThermodynamicSeaIceModel{Grid, Tim, Clk, Clo, State, Cp, Fu, Ocean, Tend}


#end

# I don't think we need a separate struct for the Winton model, however we do need to add upper layer heat capacity
# as part of the state (it's a prognostic variable)
# Right now I'm assuming the grid is a 2D grid, flat in the z-axis since the ice model "sits" on top of the ocean
# model. Thus the upper and lower temperatures are separate variables in the state. This could also be represented
# as a grid with two vertical layers and one temperature T field - this would enable generalizing to more layers later,
# but then C, hₛ, and hᵢ would have to be modified since they only have one value per horizontal cell.
function WintonThermodynamicSeaIceModel(; grid,
                                          closure = default_closure(grid),
                                          ice_density = 905.0, # kg / m³
                                          snow_density = 330.0, # kg / m³
                                          ice_heat_capacity = 2090.0 / reference_density,
                                          water_heat_capacity = 3991.0 / reference_density,
                                          latent_heat_freezing = 334000.0, # J / kg
                                          ice_thermal_conductivity = 2.03, # W / m °C
                                          snow_thermal_conductivity = 0.31, # W / m °C
                                          freezing_temp_to_salinity = 0.054, # °C per mil
                                          sea_ice_salinity = 1.0, # parts per mil
                                          freezing_temp_seawater = -1.8, # °C
                                          fusion_enthalpy = 3.3e5 / reference_density,
                                          atmosphere_temperature = -10, # ᵒC
                                          ocean_temperature = 0) # ᵒC

    # Build temperature field
    top_T_bc = ValueBoundaryCondition(atmosphere_temperature)
    bottom_T_bc = ValueBoundaryCondition(ocean_temperature)
    T_location = (Center, Center, Center)
    T_bcs_up = FieldBoundaryConditions(grid, T_location, top=top_T_bc)
    T_bcs_low = FieldBoundaryConditions(grid, T_location, bottom=bottom_T_bc)
    upper_heat_capacity = CenterField(grid)
    upper_temperature = CenterField(grid, boundary_conditions=T_bcs_up)
    lower_temperature = CenterField(grid, boundary_conditions=T_bcs_low)
    enthalpy = CenterField(grid)
    snow_thickness = CenterField(grid)
    ice_thickness = CenterField(grid)
    state = (C=upper_heat_capacity, T₁=upper_temperature, T₂=lower_temperature, H=enthalpy, hₛ=snow_thickness, hᵢ=ice_thickness)

    tendencies = (; H=CenterField(grid))
    clock = Clock{eltype(grid)}(0, 0, 1)

    return ThermodynamicSeaIceModel(grid,
                                    nothing,
                                    clock,
                                    closure,
                                    state,
                                    ice_density,
                                    snow_density,
                                    ice_heat_capacity,
                                    water_heat_capacity,
                                    latent_heat_freezing,
                                    ice_thermal_conductivity,
                                    snow_thermal_conductivity,
                                    freezing_temp_to_salinity,
                                    sea_ice_salinity,
                                    freezing_temp_seawater,
                                    fusion_enthalpy,
                                    ocean_temperature,
                                    tendencies)
end

# I think set! and the utlities can both be used as-is in ClimaSeaIce.jl. The main change will be the physics,
# and possible some time-stepping functions.


#####
##### Time-stepping
#####




#####
##### Physics
#####

end