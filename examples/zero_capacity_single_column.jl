using ClimaSeaIce
using Oceananigans
using GLMakie

# Use celsius, so
const T₀ = 273.15

function update_surface_temperature!(T_ice, h_ice, T_ocean, Q_atm, σ, k_ice, Δt)
    Nx, Ny, Nz = size(Ti)

    @inbounds begin for i = 1:Nx, j = 1:Ny
        T_ice_kelvin = T_ice[i, j, 1] + T₀

        Q_conduction = k_ice / h_ice[i, j, 1] * (T_ocean - T_ice)
        Q_emission = σ * T_ice_kelvin^4

        denominator = k_ice / h_ice[i, j, 1] - 4 * σ * T_ice_kelvin^3

        # Compute temperature increment
        ΔT = (Q_atm + Q_emission + Q_conduction) / denominator

        # Update surface temperature
        T_ice[i, j, 1] = T_ice[i, j, 1] + ΔT
    end
end

# Create a single column grid with no grid points in the vertical,
# appropriate for a zero heat capacity model.
Nx = Ny = 1
grid = RectilinearGrid(size=(Nx, Ny), x=(0, 1), y=(0, 1), topology=(Periodic, Periodic, Flat))

# Forcing and parameters
ice_conductivity = k_ice = 1e-4
atmospheric_heat_flux = Q_atm = 10 # W m⁻²
stefan_boltzman_constant = σ = 1
T_ocean = 0 # deg C

# Constant ice thickness of 1 meter
h_ice = CenterField(grid)
set!(h_ice, 1)

# Initial ice temperature of -1 Celsius
T_ice = CenterField(grid)
set!(T_ice, -1)

Δt = 1minute
Nt = 60

# Integrate forward for 1 hour
T_ice_saved = []
for n = 1:Nt
    update_surface_temperature!(T_ice, h_ice, T_ocean, Q_atm, σ, k_ice, Δt)
    push!(T_ice_saved, T_ice[1, 1, 1])
end

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, T_ice_saved)
display(fig)

