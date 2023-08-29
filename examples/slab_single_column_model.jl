using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using GLMakie

# Create a single column grid with no grid points in the vertical,
# appropriate for a zero heat capacity model.
Nx = Ny = 1
grid = RectilinearGrid(size = (Nx, Ny),
                       x = (0, 1),
                       y = (0, 1),
                       topology = (Periodic, Periodic, Flat))

extrinsic_heat_flux = Field{Center, Center, Nothing}(grid)
top_surface = MeltingRadiatingSurface(; extrinsic_heat_flux)

bottom_surface = PrescribedTemperatureSurface(0)

model = SlabSeaIceModel(grid;
                        top_surface,
                        bottom_surface,
                        ice_salinity = ConstantBulkSalinity(7), # psu
                        ice_conductivity=2.0)

simulation = Simulation(model, Δt=1minute, stop_iteration=500)

h = model.ice_thickness
set!(h, 1)

# Set initial temperature so that in and out radiation balance
ω_day = 2π / day
Qsw(t) = -300 - 50 * sin(ω_day * t)

σ = model.top_surface.stefan_boltzmann_constant
Q₀ = -300 # W m² = σ * T^4
Tu₀ = (-Q₀ / σ)^(1/4) - 273.15 # Celsius
Tu = model.top_ice_temperature
set!(Tu, Tu₀)

# Integrate forward for 1 hour
Tt = Float64[]
ht = Float64[]
times = Float64[]

function save_ice_state(sim)
    T = model.top_ice_temperature
    h = model.ice_thickness
    push!(ht, first(h))
    push!(Tt, first(T))
    push!(times, time(sim))
end

simulation.callbacks[:save] = Callback(save_ice_state)

function set_surface_shortwave!(sim)
    t = time(sim)
    Qex = sim.model.top_surface.extrinsic_heat_flux
    set!(Qex, Qsw(t))
    return nothing
end

simulation.callbacks[:fluxes] = Callback(set_surface_shortwave!)

run!(simulation)

fig = Figure(resolution=(1200, 600))
axT = Axis(fig[1, 1], title="T")
axh = Axis(fig[2, 1], title="h")
scatterlines!(axT, times ./ hour, Tt)
scatterlines!(axh, times ./ hour, ht)
display(fig)

