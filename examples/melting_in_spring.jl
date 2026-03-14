# # Melting in spring example
#
# This example simulates the melting of relatively thick sea ice in spring when
# the sun is shining. We compare bare ice with snow-covered ice to show how
# snow insulation affects melt rates.
#
# This example demonstrates how to:
#
#   * set up a one-dimensional model with multiple grid cells,
#   * prescribe spatially varying solar insolation,
#   * use `FluxFunction` for parameterized heat fluxes,
#   * add a snow layer with `snow_thermodynamics`.
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
# We generate a one-dimensional grid with 4 grid cells to model different ice
# columns subject to different solar insolation:

using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: RadiativeEmission
using CairoMakie

grid = RectilinearGrid(size=4, x=(0, 1), topology=(Periodic, Flat, Flat))

# ## Top boundary conditions
#
# We prescribe different solar insolation values for each grid cell, ranging
# from -600 to -1200 W m⁻² (negative values indicate downward/into the ice):

solar_insolation = [-600, -800, -1000, -1200] # W m⁻²
solar_insolation = reshape(solar_insolation, (4, 1, 1))

# The ice also emits longwave radiation from its top surface:

outgoing_radiation = RadiativeEmission()

# ## Sensible heat flux parameterization
#
# The sensible heat flux from the atmosphere is represented by a `FluxFunction`.
# We define the parameters for the bulk formula:

parameters = (
    transfer_coefficient     = 1e-3,  # unitless
    atmosphere_density       = 1.225, # kg m⁻³
    atmosphere_heat_capacity = 1004,  # J kg⁻¹ K⁻¹
    atmosphere_temperature   = -5,    # °C
    atmosphere_wind_speed    = 5      # m s⁻¹
)

# The flux is positive (cooling by fluxing heat upward away from the upper surface)
# when the atmosphere temperature is less than the surface temperature:

@inline function sensible_heat_flux(i, j, grid, Tu, clock, fields, parameters)
    Cs = parameters.transfer_coefficient
    ρa = parameters.atmosphere_density
    ca = parameters.atmosphere_heat_capacity
    Ta = parameters.atmosphere_temperature
    ua = parameters.atmosphere_wind_speed
    ℵ  = fields.ℵ[i, j, 1]

    return Cs * ρa * ca * ua * (Tu - Ta) * ℵ
end

aerodynamic_flux = FluxFunction(sensible_heat_flux; parameters)

# We combine all top heat fluxes into a tuple:

top_heat_flux = (outgoing_radiation, solar_insolation, aerodynamic_flux)

# ## Building the bare-ice model

model_bare = SeaIceModel(grid;
                         ice_consolidation_thickness = 0.05,
                         top_heat_flux)

set!(model_bare, h=1, ℵ=1)

# ## Building the snow-covered model
#
# We add a 20 cm layer of snow on top of the ice. Snow has a much lower
# thermal conductivity (~0.31 W/m/K) than ice (~2 W/m/K), so it acts as an
# insulating blanket that reduces the conductive flux through the slab.

snow_thermodynamics = SlabSnowThermodynamics(grid)

model_snow = SeaIceModel(grid;
                         ice_consolidation_thickness = 0.05,
                         top_heat_flux,
                         snow_thermodynamics)

set!(model_snow, h=1, ℵ=1, hs=0.2) # 20 cm of snow, no precipitation

# ## Running both simulations
#
# We run both for 30 days with a 10-minute time step, collecting time series
# for all four columns:

Δt = 10minute

# Bare-ice simulation:

simulation_bare = Simulation(model_bare, Δt=Δt, stop_time=30days)

series_bare = []

function accumulate_bare(sim)
    T = sim.model.ice_thermodynamics.top_surface_temperature
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    push!(series_bare, (time(sim),
                        h[1, 1, 1], ℵ[1, 1, 1], T[1, 1, 1],
                        h[2, 1, 1], ℵ[2, 1, 1], T[2, 1, 1],
                        h[3, 1, 1], ℵ[3, 1, 1], T[3, 1, 1],
                        h[4, 1, 1], ℵ[4, 1, 1], T[4, 1, 1]))
end

simulation_bare.callbacks[:save] = Callback(accumulate_bare)
run!(simulation_bare)

# Snow-covered simulation:

simulation_snow = Simulation(model_snow, Δt=Δt, stop_time=30days)

series_snow = []

function accumulate_snow(sim)
    m  = sim.model
    Tu = m.snow_thermodynamics.top_surface_temperature
    h  = m.ice_thickness
    ℵ  = m.ice_concentration
    hs = m.snow_thickness
    push!(series_snow, (time(sim),
                        h[1, 1, 1], ℵ[1, 1, 1], Tu[1, 1, 1], hs[1, 1, 1],
                        h[2, 1, 1], ℵ[2, 1, 1], Tu[2, 1, 1], hs[2, 1, 1],
                        h[3, 1, 1], ℵ[3, 1, 1], Tu[3, 1, 1], hs[3, 1, 1],
                        h[4, 1, 1], ℵ[4, 1, 1], Tu[4, 1, 1], hs[4, 1, 1]))
end

simulation_snow.callbacks[:save] = Callback(accumulate_snow)
run!(simulation_snow)

# ## Extracting the time series

t_bare = [d[1]  for d in series_bare]
h_bare = [[d[3*(c-1)+2] for d in series_bare] for c in 1:4]
ℵ_bare = [[d[3*(c-1)+3] for d in series_bare] for c in 1:4]
T_bare = [[d[3*(c-1)+4] for d in series_bare] for c in 1:4]

t_snow = [d[1]  for d in series_snow]
h_snow  = [[d[4*(c-1)+2] for d in series_snow] for c in 1:4]
ℵ_snow  = [[d[4*(c-1)+3] for d in series_snow] for c in 1:4]
T_snow  = [[d[4*(c-1)+4] for d in series_snow] for c in 1:4]
hs_snow = [[d[4*(c-1)+5] for d in series_snow] for c in 1:4]

# ## Visualizing the results

set_theme!(Theme(fontsize=18, linewidth=3))

fig = Figure(size=(1200, 1000))

colors = Makie.wong_colors()
labels = ["-600", "-800", "-1000", "-1200"]

# Surface temperature
axT = Axis(fig[1, 1], ylabel="Surface temperature (°C)",
           title="Surface temperature: bare (solid) vs snow-covered (dashed)")
for c in 1:4
    lines!(axT, t_bare ./ day, T_bare[c], color=colors[c], label=labels[c] * " W/m²")
    lines!(axT, t_snow ./ day, T_snow[c], color=colors[c], linestyle=:dash)
end
axislegend(axT, position=:rt)

# Ice concentration
axℵ = Axis(fig[2, 1], ylabel="Ice concentration (-)",
           title="Ice concentration: bare (solid) vs snow-covered (dashed)")
for c in 1:4
    lines!(axℵ, t_bare ./ day, ℵ_bare[c], color=colors[c])
    lines!(axℵ, t_snow ./ day, ℵ_snow[c], color=colors[c], linestyle=:dash)
end

# Ice thickness
axh = Axis(fig[3, 1], ylabel="Ice thickness (m)",
           title="Ice thickness: bare (solid) vs snow-covered (dashed)")
for c in 1:4
    lines!(axh, t_bare ./ day, h_bare[c], color=colors[c])
    lines!(axh, t_snow ./ day, h_snow[c], color=colors[c], linestyle=:dash)
end

# Snow thickness
axs = Axis(fig[4, 1], xlabel="Time (days)", ylabel="Snow thickness (m)",
           title="Snow thickness evolution")
for c in 1:4
    lines!(axs, t_snow ./ day, hs_snow[c], color=colors[c])
end

save("melting_in_spring.png", fig)
nothing # hide

# ![](melting_in_spring.png)
#
# Key observations:
#
#   * **Snow insulates:** snow-covered ice melts more slowly because the low
#     conductivity snow layer reduces the conductive flux.
#   * **Snow melts first:** the snow layer disappears before the ice starts
#     melting significantly from the top.
#   * **Warmer surface:** the snow surface is warmer than bare ice because
#     the same heat flux produces a larger temperature drop across the
#     more insulating snow+ice slab.
