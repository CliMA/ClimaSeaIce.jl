using Oceananigans
using Oceananigans.Utils: prettysummary
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using ClimaSeaIce
using ClimaSeaIce: melting_temperature
using SeawaterPolynomials: TEOS10EquationOfState
using ClimaSeaIce.ThermalBoundaryConditions: RadiativeEmission, IceWaterThermalEquilibrium
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
                       u = FieldBoundaryConditions(top=FluxBoundaryCondition(-1e-6)),
                       S = FieldBoundaryConditions(top=FluxBoundaryCondition(Qˢ)))

equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

ocean_model = HydrostaticFreeSurfaceModel(grid = ocean_grid;
                                          buoyancy, boundary_conditions, closure,
                                          tracers = (:T, :S, :e))

ocean_simulation = Simulation(ocean_model; Δt=1hour)

Nz = size(ocean_grid, 3)
So = ocean_model.tracers.S
ocean_surface_salinity = Field(So, indices=(:, :, Nz))
bottom_bc = IceWaterThermalEquilibrium(ocean_surface_salinity)

ice_model = SlabSeaIceModel(ice_grid;
                            ice_consolidation_thickness = 0.5,
                            ice_salinity = 0,
                            internal_thermal_flux = ConductiveFlux(conductivity=100),
                            top_thermal_flux = (solar_insolation, radiative_emission),
                            bottom_thermal_boundary_condition = bottom_bc,
                            bottom_thermal_flux = ice_ocean_flux)

ice_simulation = Simulation(ice_model, Δt=1hour)


using SeawaterPolynomials: thermal_expansion, haline_contraction

# Double stratification
N²θ = 0
T₀ = 0
S₀ = 30
g = ocean_model.buoyancy.model.gravitational_acceleration
α = thermal_expansion(T₀, S₀, 0, equation_of_state)
dTdz = N²θ / (α * g)

N²S = 1e-4
β = haline_contraction(T₀, S₀, 0, equation_of_state)
dSdz = N²S / (β * g)

Tᵢ(x, y, z) = T₀ + z * dTdz
Sᵢ(x, y, z) = S₀ - z * dSdz

set!(ocean_model, S=Sᵢ, T=T₀)

coupled_model = IceOceanModel(ice_simulation, ocean_simulation)

t = Float64[]
hi = Float64[]
αi = Float64[]
Tst = Float64[]

To = []
So = []
eo = []

Nz = size(ocean_grid, 3)

while(time(ocean_simulation) < 100days)
    time_step!(coupled_model, 20minutes)

    if mod(iteration(ocean_simulation), 10) == 0
        push!(t, time(ocean_simulation))

        h = ice_model.ice_thickness
        α = ice_model.ice_concentration
        Ts = ice_model.top_surface_temperature
        push!(hi, first(h))
        push!(αi, first(α))
        push!(Tst, first(Ts))

        T = ocean_model.tracers.T
        S = ocean_model.tracers.S
        e = ocean_model.tracers.e

        push!(To, deepcopy(interior(T, 1, 1, :)))
        push!(So, deepcopy(interior(S, 1, 1, :)))
        push!(eo, deepcopy(interior(e, 1, 1, :)))

        @info string("Iter: ", iteration(ocean_simulation),
                     ", t: ", prettytime(ocean_simulation),
                     ", Tₒ: ", prettysummary(ocean_model.tracers.T[1, 1, Nz]),
                     ", Sₒ: ", prettysummary(ocean_model.tracers.S[1, 1, Nz]),
                     ", h: ", prettysummary(hi[end]),
                     ", α: ", prettysummary(αi[end]),
                     ", Tᵢ: ", prettysummary(Tst[end]),
                     ", α h: ", prettysummary(αi[end] * hi[end]))
    end
end

Ts = map(T -> T[Nz], To)
Ss = map(T -> T[Nz], So)

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(2400, 1200))

axhi = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")
axαi = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice concentration")
axTo = Axis(fig[3, 1], xlabel="Time (days)", ylabel="Surface temperatures (ᵒC)")
axSo = Axis(fig[4, 1], xlabel="Time (days)", ylabel="Ocean surface salinity (psu)")

lines!(axhi, t ./ day, hi)
lines!(axαi, t ./ day, αi)
lines!(axTo, t ./ day, Ts, label="Ocean")
lines!(axTo, t ./ day, Tst, label="Ice")
axislegend(axTo)

lines!(axSo, t ./ day, Ss)

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

vlines!(axhi, tn, color=(:yellow, 0.6))
vlines!(axαi, tn, color=(:yellow, 0.6))
vlines!(axTo, tn, color=(:yellow, 0.6))
vlines!(axSo, tn, color=(:yellow, 0.6))

lines!(axTz, Tn, z)
lines!(axSz, Sn, z)

xlims!(axTz, -2, 6)
# xlims!(axSz, 26, 31)

colsize!(fig.layout, 2, Relative(0.2))
colsize!(fig.layout, 3, Relative(0.2))

display(fig)

record(fig, "freezing_and_melting.mp4", 1:Nt, framerate=12) do nn
    n[] = nn
end

