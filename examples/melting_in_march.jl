using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using GLMakie

# Generate a 0D grid for a single column slab model 
grid = RectilinearGrid(size=4, x=(0, 1), topology=(Periodic, Flat, Flat))

# Build a model of an ice slab that has internal conductive fluxes
# and that emits radiation from its top surface.
solar_insolation = [-600, -800, -1000, -1200] # W m⁻²
solar_insolation = reshape(solar_insolation, (4, 1, 1))
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

top_heat_flux = (outgoing_radiation, solar_insolation, aerodynamic_flux)
model = SlabSeaIceModel(grid;
                        ice_consolidation_thickness = 0.01, # m
                        top_heat_flux)
set!(model, h=1)

simulation = Simulation(model, Δt=10minute, stop_time=30days)

# Accumulate data
timeseries = []

function accumulate_timeseries(sim)
    T = model.top_surface_temperature
    h = model.ice_thickness
    push!(timeseries, (time(sim),
                       h[1, 1, 1], T[1, 1, 1],
                       h[2, 1, 1], T[2, 1, 1],
                       h[3, 1, 1], T[3, 1, 1],
                       h[4, 1, 1], T[4, 1, 1]))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

run!(simulation)

# Extract and visualize data

t = [datum[1] for datum in timeseries]
h1 = [datum[2] for datum in timeseries]
T1 = [datum[3] for datum in timeseries]
h2 = [datum[4] for datum in timeseries]
T2 = [datum[5] for datum in timeseries]
h3 = [datum[6] for datum in timeseries]
T3 = [datum[7] for datum in timeseries]
h4 = [datum[8] for datum in timeseries]
T4 = [datum[9] for datum in timeseries]

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(1000, 800))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top temperature (ᵒC)")
axh = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")

lines!(axT, t / day, T1)
lines!(axh, t / day, h1)

lines!(axT, t / day, T2)
lines!(axh, t / day, h2)

lines!(axT, t / day, T3)
lines!(axh, t / day, h3)

lines!(axT, t / day, T4)
lines!(axh, t / day, h4)

display(fig)

