# # Freezing of a lake example
#
# In this example we simulate the freezing of a lake in winter. The lake is
# represented by four grid points that start at 1°C and are cooled down by an
# atmosphere with temperatures of -20, -10, -5, and 0°C. The lake is 10 m deep
# and not subjected to radiative transfer (this lake is covered in a wind tunnel
# under which we blow some cold air). This example demonstrates
#
#   * How to couple a simple lake model with sea ice thermodynamics.
#   * How to use `FluxFunction` for complex boundary conditions.
#   * How to track energy budgets.
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.
#
# ```julia
# using Pkg
# pkg"add Oceananigans, ClimaSeaIce, CairoMakie"
# ```
#
# ## The physical domain
#
# We generate a one-dimensional grid with 4 grid cells to model different lake
# columns subject to different atmospheric temperatures:

using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: latent_heat
using CairoMakie

grid = RectilinearGrid(size=4, x=(0, 1), topology=(Periodic, Flat, Flat))

# ## Atmospheric forcing
#
# The sensible heat flux from the atmosphere is represented by a `FluxFunction`.
# We define the atmospheric parameters:

atmosphere = (
    transfer_coefficient     = 1e-3,  # unitless
    atmosphere_density       = 1.225, # kg m⁻³
    atmosphere_heat_capacity = 1004,  # J kg⁻¹ K⁻¹
    atmosphere_temperature   = [-20, -10, -5, 0], # °C
    atmosphere_wind_speed    = 5,     # m s⁻¹
    atmosphere_ice_flux      = [0.0, 0.0, 0.0, 0.0], # W m⁻²
)
    
# The flux is positive (cooling by fluxing heat upward away from the upper surface)
# when the atmosphere temperature is less than the surface temperature:

@inline function sensible_heat_flux(i, j, grid, Tᵤ, clock, fields, atmosphere)
    Cₛ = atmosphere.transfer_coefficient
    ρₐ = atmosphere.atmosphere_density
    cₐ = atmosphere.atmosphere_heat_capacity
    Tₐ = atmosphere.atmosphere_temperature[i]
    uₐ = atmosphere.atmosphere_wind_speed
    ℵ  = fields.ℵ[i, j, 1]
    Qₐ = atmosphere.atmosphere_ice_flux

    Qₐ[i] = ifelse(ℵ == 0, zero(grid), Cₛ * ρₐ * cₐ * uₐ * (Tᵤ - Tₐ))

    return Qₐ[i] 
end

# ## Lake model
#
# We also evolve a simple bucket freshwater lake that cools down and freezes from
# above, generating fresh sea ice (or lake ice in this case). We set the time step
# of the lake to 10 minutes, which will also be used for the sea ice model:

lake = (
    lake_density         = 1000, # kg m⁻³
    lake_heat_capacity   = 4000,  # J kg⁻¹ K⁻¹
    lake_temperature     = [1.0, 1.0, 1.0, 1.0], # °C
    lake_depth           = 10, # m
    lake_ice_flux        = [0.0, 0.0, 0.0, 0.0], # W m⁻²
    atmosphere_lake_flux = [0.0, 0.0, 0.0, 0.0], # W m⁻²
    Δt                   = 10minutes
)

# ## Bottom heat flux function
#
# We build a flux function that serves three purposes:
# 1. Computing the cooling of the lake caused by the atmosphere
# 2. If the lake temperature is low enough, freezing the lake from above
# 3. Adding the heat flux to the bottom of the ice

@inline function advance_lake_and_frazil_flux(i, j, grid, Tuᵢ, clock, fields, parameters)
    atmos = parameters.atmosphere
    lake  = parameters.lake

    Cₛ = atmos.transfer_coefficient
    ρₐ = atmos.atmosphere_density
    cₐ = atmos.atmosphere_heat_capacity
    Tₐ = atmos.atmosphere_temperature[i]
    uₐ = atmos.atmosphere_wind_speed
    Tₒ = lake.lake_temperature
    cₒ = lake.lake_heat_capacity
    ρₒ = lake.lake_density
    Δ  = lake.lake_depth
    Δt = lake.Δt
    ℵ  = fields.ℵ[i, j, 1]
    
    Qᵗ = lake.atmosphere_lake_flux
    Qᵇ = lake.lake_ice_flux
    
    Qₐ = Cₛ * ρₐ * cₐ * uₐ * (Tₐ - Tₒ[i]) * (1 - ℵ)
    Tⁿ = Tₒ + Qₐ / (ρₒ * cₒ) * Δt
    Qᵢ = ρₒ * cₒ * (Tⁿ - 0) / Δt * Δ # W m⁻²
    Qᵢ = min(Qᵢ, zero(Qᵢ))

    Tₒ[i] = ifelse(Qᵢ == 0, Tⁿ, zero(Qᵢ))
    Qᵗ[i] = Qₐ
    Qᵇ[i] = Qᵢ

    return Qᵢ
end

# ## Building the model
#
# We assemble the model with the top and bottom heat flux functions:

top_heat_flux    = FluxFunction(sensible_heat_flux; parameters=atmosphere)
bottom_heat_flux = FluxFunction(advance_lake_and_frazil_flux; parameters=(; lake, atmosphere))

model = SeaIceModel(grid;
                    ice_consolidation_thickness = 0.05, # m
                    top_heat_flux, 
                    bottom_heat_flux)

# We initialize all columns with open water (0 thickness and 0 concentration):

set!(model, h=0, ℵ=0)

# ## Running a simulation
#
# We set up a simulation that runs for 20 days:

simulation = Simulation(model, Δt=lake.Δt, stop_time=20days)

# ## Collecting data
#
# We set up callbacks to accumulate time series data and energy budgets:

timeseries = []

function accumulate_timeseries(sim)
    T  = model.ice_thermodynamics.top_surface_temperature
    h  = model.ice_thickness
    ℵ  = model.ice_concentration
    To = lake.lake_temperature
    push!(timeseries, (time(sim),
                       h[1, 1, 1], ℵ[1, 1, 1], T[1, 1, 1],
                       h[2, 1, 1], ℵ[2, 1, 1], T[2, 1, 1],
                       h[3, 1, 1], ℵ[3, 1, 1], T[3, 1, 1],
                       h[4, 1, 1], ℵ[4, 1, 1], T[4, 1, 1],
                       To[1], To[2], To[3], To[4]))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

# We also accumulate energy for budget analysis:

Ei = []
Qa = []
Ql = []

function accumulate_energy(sim)
    h  = sim.model.ice_thickness
    ℵ  = sim.model.ice_concentration
    PT = sim.model.ice_thermodynamics.phase_transitions
    ℰ  = latent_heat(PT, 0) # ice is at 0°C

    push!(Ei, deepcopy(@. - h * ℵ * ℰ))
    push!(Qa, deepcopy(atmosphere.atmosphere_ice_flux))
    push!(Ql, deepcopy(lake.lake_ice_flux))
end

simulation.callbacks[:energy] = Callback(accumulate_energy)

run!(simulation)

# ## Visualizing the results
#
# We extract the time series data for all four columns:

t  = [datum[1]  for datum in timeseries]
h1 = [datum[2]  for datum in timeseries]
ℵ1 = [datum[3]  for datum in timeseries]
T1 = [datum[4]  for datum in timeseries]
h2 = [datum[5]  for datum in timeseries]
ℵ2 = [datum[6]  for datum in timeseries]
T2 = [datum[7]  for datum in timeseries]
h3 = [datum[8]  for datum in timeseries]
ℵ3 = [datum[9]  for datum in timeseries]
T3 = [datum[10] for datum in timeseries]
h4 = [datum[11] for datum in timeseries]
ℵ4 = [datum[12] for datum in timeseries]
T4 = [datum[13] for datum in timeseries]
L1 = [datum[14] for datum in timeseries]
L2 = [datum[15] for datum in timeseries]
L3 = [datum[16] for datum in timeseries]
L4 = [datum[17] for datum in timeseries]

# And visualize the evolution of ice thickness, concentration, surface temperature,
# and lake temperature:

set_theme!(Theme(fontsize=18, linewidth=3))

fig = Figure(size=(1000, 900))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top temperature (ᵒC)")
axℵ = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice concentration (-)")
axh = Axis(fig[3, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")
axL = Axis(fig[4, 1], xlabel="Time (days)", ylabel="Lake temperature (ᵒC)")

lines!(axT, t ./ day, T1)
lines!(axℵ, t ./ day, ℵ1)
lines!(axh, t ./ day, h1)
lines!(axL, t ./ day, L1)

lines!(axT, t ./ day, T2)
lines!(axℵ, t ./ day, ℵ2)
lines!(axh, t ./ day, h2)
lines!(axL, t ./ day, L2)

lines!(axT, t ./ day, T3)
lines!(axℵ, t ./ day, ℵ3)
lines!(axh, t ./ day, h3)
lines!(axL, t ./ day, L3)

lines!(axT, t ./ day, T4)
lines!(axℵ, t ./ day, ℵ4)
lines!(axh, t ./ day, h4)
lines!(axL, t ./ day, L4)

save("freezing_in_winter.png", fig)
nothing # hide

# ![](freezing_in_winter.png)

# ## Energy budget analysis
#
# We also visualize the energy budget to verify conservation:

Ei1 = [datum[1]  for datum in Ei]
Qa1 = [datum[1]  for datum in Qa]
Ql1 = [datum[1]  for datum in Ql]
Ei2 = [datum[2]  for datum in Ei]
Qa2 = [datum[2]  for datum in Qa]
Ql2 = [datum[2]  for datum in Ql]
Ei3 = [datum[3]  for datum in Ei]
Qa3 = [datum[3]  for datum in Qa]
Ql3 = [datum[3]  for datum in Ql]
Ei4 = [datum[4]  for datum in Ei]
Qa4 = [datum[4]  for datum in Qa]
Ql4 = [datum[4]  for datum in Ql]

fig = Figure(size=(1000, 900))

axE = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Sea ice energy (J)")
axA = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Atmosphere heat flux (W m⁻²)")
axL = Axis(fig[3, 1], xlabel="Time (days)", ylabel="Lake heat flux (W m⁻²)")
axB = Axis(fig[4, 1], xlabel="Time (days)", ylabel="Heat budget residual (W m⁻²)")

dEi1 = (Ei1[2:end] - Ei1[1:end-1]) ./ 10minutes
dEi2 = (Ei2[2:end] - Ei2[1:end-1]) ./ 10minutes
dEi3 = (Ei3[2:end] - Ei3[1:end-1]) ./ 10minutes
dEi4 = (Ei4[2:end] - Ei4[1:end-1]) ./ 10minutes
tpE  = t[2:end] 

lines!(axE, t ./ day, Ei1)
lines!(axA, t ./ day, Qa1)
lines!(axL, t ./ day, Ql1)

lines!(axE, t ./ day, Ei2)
lines!(axA, t ./ day, Qa2)
lines!(axL, t ./ day, Ql2)

lines!(axE, t ./ day, Ei3)
lines!(axA, t ./ day, Qa3)
lines!(axL, t ./ day, Ql3)

lines!(axE, t ./ day, Ei4)
lines!(axA, t ./ day, Qa4)
lines!(axL, t ./ day, Ql4)

# The budget residual should be close to zero
lines!(axB, tpE ./ day, dEi1 .- (.- Qa1 .+ Ql1)[2:end])
lines!(axB, tpE ./ day, dEi2 .- (.- Qa2 .+ Ql2)[2:end])
lines!(axB, tpE ./ day, dEi3 .- (.- Qa3 .+ Ql3)[2:end])

save("energy_budget.png", fig)
nothing # hide

# ![](energy_budget.png)
