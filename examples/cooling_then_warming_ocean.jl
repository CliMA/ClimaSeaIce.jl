using Oceananigans
using Oceananigans.Utils: prettysummary
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using ClimaSeaIce
using ClimaSeaIce: melting_temperature
using SeawaterPolynomials: TEOS10EquationOfState
using ClimaSeaIce.ThermalBoundaryConditions: RadiativeEmission
using GLMakie

import Oceananigans.Simulations: time_step!, time

include("ice_ocean_model.jl")

Δt = 4hours
mixed_layer_depth = hₒ = 100
single_column_ocean = true

ice_grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

if single_column_ocean
    ocean_grid = RectilinearGrid(size=20, z=(-hₒ, 0), topology=(Flat, Flat, Bounded))
    closure = CATKEVerticalDiffusivity()
else
    ocean_grid = RectilinearGrid(size=1, z=(-hₒ, 0), topology=(Flat, Flat, Bounded))
    closure = nothing
end

# Top boundary conditions:
#   - outgoing radiative fluxes emitted from surface
#   - incoming shortwave radiation starting after 40 days

radiative_emission = RadiativeEmission()
top_ocean_heat_flux = Qᵀ = Field{Center, Center, Nothing}(ocean_grid)
ice_ocean_flux = Field{Center, Center, Nothing}(ice_grid)

solar_insolation = I₀ = Field{Center, Center, Nothing}(ocean_grid)
compute_solar_insolation!(sim) = (time(sim) > 40days) && (@inbounds I₀[1, 1, 1] = -600) # W m⁻²

# Generate a zero-dimensional grid for a single column slab model 

top_salt_flux = Qˢ = Field{Center, Center, Nothing}(ocean_grid)
boundary_conditions = (T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                       S = FieldBoundaryConditions(top=FluxBoundaryCondition(Qˢ)))

equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

ocean_model = HydrostaticFreeSurfaceModel(grid = ocean_grid;
                                          buoyancy, boundary_conditions, closure,
                                          tracers = (:T, :S, :e))

ocean_simulation = Simulation(ocean_model; Δt=1hour)


ice_model = SlabSeaIceModel(ice_grid;
                            minimum_ice_thickness = 0.01,
                            top_thermal_flux = (solar_insolation, radiative_emission),
                            bottom_thermal_flux = ice_ocean_flux)

ice_simulation = Simulation(ice_model, Δt=1hour)


using SeawaterPolynomials: thermal_expansion, haline_contraction

# Double stratification
N²θ = 0
T₀ = -1
S₀ = 30
g = ocean_model.buoyancy.model.gravitational_acceleration
α = thermal_expansion(T₀, S₀, 0, equation_of_state)
dTdz = N²θ / (α * g)

N²S = 1e-4
β = haline_contraction(T₀, S₀, 0, equation_of_state)
dSdz = N²S / (β * g)

Tᵢ(x, y, z) = T₀ + z * dTdz
Sᵢ(x, y, z) = S₀ - z * dSdz

set!(ocean_model, S=Sᵢ, T=-1.6) #Tᵢ)

coupled_model = IceOceanModel(ice_simulation, ocean_simulation)

t = Float64[]
hi = Float64[]
αi = Float64[]

To = []
So = []
eo = []

Nz = size(ocean_grid, 3)

for i = 1:1000
    time_step!(coupled_model, 20minutes)

    if mod(i, 1) == 0
        push!(t, time(ocean_simulation))

        h = ice_model.ice_thickness
        α = ice_model.ice_concentration
        push!(hi, first(h))
        push!(αi, first(α))

        T, S, e = ocean_model.tracers
        push!(To, deepcopy(interior(T, 1, 1, :)))
        push!(So, deepcopy(interior(S, 1, 1, :)))
        push!(eo, deepcopy(interior(e, 1, 1, :)))

        @info string("Iter: ", iteration(ocean_simulation),
                     ", t: ", prettytime(ocean_simulation),
                     ", Tₒ: ", prettysummary(ocean_model.tracers.T[1, 1, Nz]),
                     ", Sₒ: ", prettysummary(ocean_model.tracers.S[1, 1, Nz]),
                     ", h: ", prettysummary(hi[end]),
                     ", α: ", prettysummary(αi[end]),
                     ", α h: ", prettysummary(αi[end] * hi[end]))
    end
end

Ts = map(T -> T[Nz], To)
Ss = map(T -> T[Nz], So)

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(2400, 1200))

axh = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")
axα = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice concentration")
axT = Axis(fig[3, 1], xlabel="Time (days)", ylabel="Ocean surface temperature (ᵒC)")
axS = Axis(fig[4, 1], xlabel="Time (days)", ylabel="Ocean surface salinity (psu)")

lines!(axh, t ./ day, hi)
lines!(axα, t ./ day, αi)
lines!(axT, t ./ day, Ts)
lines!(axS, t ./ day, Ss)

axTz = Axis(fig[1:4, 2], xlabel="T", ylabel="z (m)")
axSz = Axis(fig[1:4, 3], xlabel="S", ylabel="z (m)")

z = znodes(ocean_model.tracers.T)

Nt = length(t)
slider = Slider(fig[5, 1:3], range=1:Nt, startvalue=1)
n = slider.value

Tn = @lift To[$n]
Sn = @lift So[$n]
en = @lift max.(1e-6, eo[$n])

tn = @lift t[$n] / day

vlines!(axh, tn, color=(:black, 0.6))
vlines!(axT, tn, color=(:black, 0.6))
vlines!(axS, tn, color=(:black, 0.6))

lines!(axTz, Tn, z)
lines!(axSz, Sn, z)

xlims!(axTz, -2, 6)
# xlims!(axSz, 26, 31)

colsize!(fig.layout, 2, Relative(0.2))
colsize!(fig.layout, 3, Relative(0.2))

display(fig)

# record(fig, "freezing_and_melting.mp4", 1:Nt, framerate=12) do nn
#     n[] = nn
# end

