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

m = 0.054 #ᵒ psu⁻¹
@show Tm = - m * 35 # ᵒC
bottom_surface = PrescribedTemperatureSurface(Tm)

model = SlabSeaIceModel(grid;
                        top_surface,
                        bottom_surface,
                        reference_fusion_enthalpy = 333e3,
                        ice_salinity = ConstantBulkSalinity(7), # psu
                        ice_conductivity = 2.0)

simulation = Simulation(model, Δt=10minute, stop_time=30days)

# Set initial ice thickness and surface temperature
h = model.ice_thickness
Tu = model.top_ice_temperature
set!(h, 2)
set!(Tu, -10)

# Parameters: latitude of the column and ice albedo
φ = 81  # degrees
α = 0.7 # albedo

# Typical normal component of clear sky solar radiation
maximum_clear_sky_radiative_flux = -1025 # W m⁻²

# Earth's obliquity or tilt
obliquity = -23.44 # degrees

# Solar decliniation angle
δ(t) = asin(sind(obliquity) * sin(π + 2π/365 * t/day))

# Solar hour angle
h(t) = 15π / 180 * (t/hour - 12)

# Solar zenith angle
cosine_zenith_angle(t) = sin(φ) * sin(δ(t)) + cos(φ) * cos(δ(t)) * sin(h(t))

# Put it all together: Incoming shortwave heat flux
Q_shortwave(t) = (1 - α) * maximum_clear_sky_radiative_flux * max(0, cosine_zenith_angle(t))

# Model for turbulent fluxes
C = 1e-3   # transfer coefficient
U = 10     # m s⁻¹ light winds, if you will
ρₐ = 1.225 # typical density of air
cₐ = 1004  # heat capacity of dry air
Tₐ = -5    # try this
Q_sensible(Tᵢ, Tₐ) = C * ρₐ * cₐ * U * (Tᵢ - Tₐ)

Qxt = Float64[]

function compute_extrinsic_heat_flux!(sim)
    t = time(sim)
    Tᵢ = first(sim.model.top_ice_temperature)
    Qx = Q_shortwave(t) + Q_sensible(Tᵢ, Tₐ)
    push!(Qxt, Qx)
    set!(sim.model.top_surface.extrinsic_heat_flux, Qx)
    return nothing
end

simulation.callbacks[:fluxes] = Callback(compute_extrinsic_heat_flux!)

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

run!(simulation)

set_theme!(Theme(fontsize=24, linewidth=4))
fig = Figure(resolution=(1200, 600))
axQ = Axis(fig[1, 1], xlabel="Time (hours)", ylabel="Shortwave heat flux (W m⁻²)")
axT = Axis(fig[2, 1], xlabel="Time (hours)", ylabel="Top surface temperature (ᵒC)")
axh = Axis(fig[3, 1], xlabel="Time (hours)", ylabel="Ice thickness (m)")
lines!(axQ, times ./ hour, Qxt)
lines!(axT, times ./ hour, Tt)
lines!(axh, times ./ hour, ht)
display(fig)

