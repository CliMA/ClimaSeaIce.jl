using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using GLMakie

# Generate a 0D grid for a single column slab model 
grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# Build a model of an ice slab that has internal conductive fluxes
# and that emits radiation from its top surface.
solar_insolation = -500 # W m⁻²
outgoing_radiation = RadiativeEmission()

# Define a FluxFunction representing a sensible heat flux
parameters = (
    transfer_coefficient     = 1e-3,  # Unitless
    atmosphere_density       = 1.225, # kg m⁻³
    atmosphere_heat_capacity = 1004,  # 
    atmosphere_temperature   = -5,    # ᵒC
    atmosphere_wind_speed    = 5      # m s⁻¹
)

@inline function sensible_heat_flux(i, j, grid, Tᵤ, clock, fields, parameters)
    Cₛ = parameters.transfer_coefficient
    ρₐ = parameters.atmosphere_density
    cₐ = parameters.atmosphere_heat_capacity
    Tₐ = parameters.atmosphere_temperature
    uₐ = parameters.atmosphere_wind_speed

    # Flux is positive (cooling by fluxing heat up away from upper surface)
    # when Tₐ < Tᵤ:
    return Cₛ * ρₐ * cₐ * uₐ * (Tᵤ - Tₐ)
end

aerodynamic_flux = FluxFunction(sensible_heat_flux; parameters)

top_thermal_flux = (outgoing_radiation, solar_insolation, aerodynamic_flux)
model = SlabSeaIceModel(grid; top_thermal_flux)
set!(model, h=1)

simulation = Simulation(model, Δt=10minute, stop_time=1day)

# Accumulate data
timeseries = []

function accumulate_timeseries(sim)
    T = model.top_temperature
    h = model.ice_thickness
    push!(timeseries, (time(sim), first(h), first(T)))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

run!(simulation)

# Extract and visualize data

t = [datum[1] for datum in timeseries]
h = [datum[2] for datum in timeseries]
T = [datum[3] for datum in timeseries]

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(1000, 800))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top temperature (ᵒC)")
axh = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")

lines!(axT, t / day, T)
lines!(axh, t / day, h)

display(fig)

