# Information about the diurnal temperature fluctuations and cooling trends at the atmosphere boundary
top_T_initial = -5
top_T_amplitude = 5
top_T_slope = -0.5/86400 #degree C/sec

# Information about ocean cooling
bottom_T_initial = 0
bottom_T_slope = -0.5/86400 #degree C/sec 

# Information about time steps
dt = 1000
tf = 864000
time = collect(0:dt:tf)

# Calculate BCs
top_T_BC=top_T_slope*time.+top_T_amplitude*sin.((2*pi*time)/86400).+top_T_initial
bottom_T_BC=bottom_T_slope*time.+bottom_T_initial

# Plot boundary conditions
using Plots
plot(time, [top_T_BC bottom_T_BC], title="Boundary Conditions", label=["Atmosphere Temp." "Ocean Temp."], linewidth=3, xlabel="Time (s)", ylabel="Temperature (C)")
