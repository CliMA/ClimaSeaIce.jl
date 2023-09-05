using Oceananigans
using Oceananigans.Utils: prettysummary
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce: melting_temperature
using SeawaterPolynomials: TEOS10EquationOfState
using ClimaSeaIce.ThermalBoundaryConditions: RadiativeEmission
using GLMakie

import Oceananigans.Simulations: time_step!, time

include("coupled_model.jl")

Δt = 4hours
mixed_layer_depth = hₒ = 100

ice_grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
ocean_grid = RectilinearGrid(size=1, z=(-hₒ, 0), topology=(Flat, Flat, Bounded))

# Top boundary conditions:
#   - outgoing radiative fluxes emitted from surface
#   - incoming shortwave radiation starting after 40 days

radiative_emission = RadiativeEmission()
top_ocean_heat_flux = Qᵀ = Field{Center, Center, Nothing}(ocean_grid)
ice_ocean_flux = Field{Center, Center, Nothing}(ice_grid)

solar_insolation = I₀ = Field{Center, Center, Nothing}(ocean_grid)
compute_solar_insolation!(sim) = (time(sim) > 40days) && (@inbounds I₀[1, 1, 1] = -400) # W m⁻²

# Generate a zero-dimensional grid for a single column slab model 

top_salt_flux = Qˢ = Field{Center, Center, Nothing}(ocean_grid)
boundary_conditions = (T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                       S = FieldBoundaryConditions(top=FluxBoundaryCondition(Qˢ)))

equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

ocean_model = HydrostaticFreeSurfaceModel(grid = ocean_grid; buoyancy, boundary_conditions,
                                          velocities = PrescribedVelocityFields(),
                                          tracers = (:T, :S))

ocean_simulation = Simulation(ocean_model; Δt=1hour)


ice_model = SlabSeaIceModel(ice_grid;
                            top_thermal_flux = (solar_insolation, radiative_emission),
                            bottom_thermal_flux = ice_ocean_flux)

ice_simulation = Simulation(ice_model, Δt=1hour)

set!(ocean_model, S=35, T=0)

coupled_model = IceOceanModel(ice_simulation, ocean_simulation)

To = Float64[]
So = Float64[]
hi = Float64[]
t = Float64[]

for i = 1:3000
    time_step!(coupled_model, 1hour)

    push!(To, first(ocean_model.tracers.T))
    push!(So, first(ocean_model.tracers.S))
    push!(hi, first(ice_model.ice_thickness))
    push!(t, time(ocean_simulation))

    @info string("Iter: ", iteration(ocean_simulation),
                 ", time: ", prettytime(ocean_simulation),
                 ", ocean temperature: ", prettysummary(first(ocean_model.tracers.T)),
                 ", ocean salinity: ", prettysummary(first(ocean_model.tracers.S)),
                 ", ice thickness: ", prettysummary(first(ice_model.ice_thickness)))
end

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(1200, 1200))

axh = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")
axT = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ocean temperature (ᵒC)")
axS = Axis(fig[3, 1], xlabel="Time (days)", ylabel="Ocean salinity (psu)")

lines!(axh, t ./ day, hi)
lines!(axT, t ./ day, To)
lines!(axS, t ./ day, So)

display(fig)

