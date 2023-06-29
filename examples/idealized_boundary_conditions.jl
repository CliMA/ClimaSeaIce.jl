using Oceananigans.Units
using GLMakie

# This file creates idealized ice-interface temperatures
# for driving an ice column model.

# Information about the diurnal temperature fluctuations and cooling trends at the atmosphere boundary
top_T_initial = -5
top_T_amplitude = 5
top_T_slope = -0.5 / day # ᵒC s⁻¹

# Information about ocean cooling
bottom_T_initial = 0
bottom_T_slope = -0.5 / day # ᵒC s⁻¹

# Information about time steps
dt = 1000
tf = 864000
time = 0:dt:tf

# Calculate BCs
top_T_BC = @. top_T_slope * time + top_T_amplitude * sin(2π*time/day) + top_T_initial
bottom_T_BC = @. bottom_T_slope * time + bottom_T_initial

# Plot boundary conditions

set_theme!(Theme(fontsize=24, linewidth=3))
fig = Figure()
ax = Axis(fig[1, 1], title="Boundary Conditions", xlabel="Time (s)", ylabel="Temperature (ᵒC)")
lines!(ax, time, top_T_BC, label="Air-ice surface temperature")
lines!(ax, time, bottom_T_BC, label="Ice-ocean temperature")
axislegend(ax)
     
display(fig)

