using Oceananigans
using Oceananigans.Units
using ClimaSeaIce
using ClimaSeaIce: melting_temperature
using SeawaterPolynomials: TEOS10EquationOfState
using ClimaSeaIce.ThermalBoundaryConditions: RadiativeEmission
using GLMakie

import Oceananigans.Simulations: time_step!

# Ocean density and heat capacity
const ρₒ = 1024
const cₒ = 3991
Δt = 4hours

# Generate a zero-dimensional grid for a single column slab model 
ocean_grid = RectilinearGrid(size=1, z=(-100, 0), topology=(Flat, Flat, Bounded))

Qᵀ = Field{Center, Center, Nothing}(ocean_grid)
Qˢ = Field{Center, Center, Nothing}(ocean_grid)
T_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ))
S_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Qˢ))

equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

ocean_model = HydrostaticFreeSurfaceModel(grid = ocean_grid; buoyancy,
                                          velocities = PrescribedVelocityFields(),
                                          tracers = (:T, :S),
                                          boundary_conditions = (; T=T_bcs, S=S_bcs))

ocean_simulation = Simulation(ocean_model; Δt=1hour)

# Build a model of an ice slab that has internal conductive fluxes
# and that emits radiation from its top surface.
ice_grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

# Ice-ocean heat flux
Qₒ = Field{Center, Center, Nothing}(ice_grid)

ice_model = SlabSeaIceModel(grid;
                            top_thermal_flux = RadiativeEmission(),
                            bottom_thermal_flux = Qₒ)

ice_simulation = Simulation(ice_model, Δt=1hour)

struct IceOceanSimulation{I, O}
    ice :: I
    ocean :: O
end

set!(ocean_model, S=35, T=0)

function compute_air_sea_flux!(coupled_sim)
    ocean_sim = coupled_sim.ocean
    ice_sim = coupled_sim.ice
    radiation = ice_sim.model.external_thermal_fluxes.top

    ocean_grid = ocean_sim.model.grid
    T = ocean_sim.model.tracers.T
    Nz = size(ocean_grid, 3)
    T₀ = interior(T, :, :, Nz)
    Qᵀ = T.boundary_conditions.top.condition

    ϵ = 1.0 # ocean emissivity
    σ = radiation.stefan_boltzmann_constant
    Tᵣ = radiation.reference_temperature

    Qᵀi = interior(Qᵀ, :, :, 1)
    @. Qᵀi = ϵ * σ * (T₀ + Tᵣ)^4 / (ρₒ * cₒ)

    h = @inbounds ice_sim.model.ice_thickness[1, 1, 1]
    @inbounds Qᵀi[1, 1, 1] = ifelse(h > 0, zero(grid), Qᵀi[1, 1, 1])

    return nothing
end

coupled_sim = IceOceanSimulation(ice_simulation, ocean_simulation)
compute_air_sea_flux!(coupled_sim)

function time_step!(coupled_sim::IceOceanSimulation, Δt)
    ocean_sim = coupled_sim.ocean
    ice_sim   = coupled_sim.ice
    liquidus  = ice_sim.model.phase_transitions.liquidus
    ice_sim.Δt = Δt
    ocean_sim.Δt = Δt

    compute_air_sea_flux!(coupled_sim)
    time_step!(ocean_sim)

    # Check ocean temperature
    hₒ = ocean_sim.model.grid.Lz # mixed layer depth
    Qₒ = ice_sim.model.external_thermal_fluxes.bottom
    Tₒ = ocean_sim.model.tracers.T
    Sₒ = ocean_sim.model.tracers.S

    @inbounds begin
        T₁ = Tₒ[1, 1, 1]
        S₁ = Sₒ[1, 1, 1]

        # Compute total latent heat (per unit area) and latent heat flux
        Tₘ = melting_temperature(liquidus, S₁)
        δE = ρₒ * cₒ * (T₁ - Tₘ) # > 0
        δQ = δE * hₒ / Δt        # > 0 (we are warming the ocean)

        # Clip 
        Qₒ[1, 1, 1] = min(zero(grid), δQ)
        Tₒ[1, 1, 1] = max(Tₘ, Tₒ[1, 1, 1])
    end

    time_step!(ice_sim)

    # TODO after ice time-step:
    #   - Adjust ocean temperature if the ice completely melts?
    #   - Adjust salinity (or compute salinity flux?)
    
    return nothing
end

To = Float64[]
hi = Float64[]
t = Float64[]

for i = 1:2000
    time_step!(coupled_sim, 1hour)

    push!(To, first(ocean_model.tracers.T))
    push!(hi, first(ice_model.ice_thickness))
    push!(t, time(ocean_simulation))

    @info string("Iter: ", iteration(ocean_simulation),
                 ", time: ", prettytime(ocean_simulation),
                 ", ocean temperature: ", first(ocean_model.tracers.T),
                 ", ice thickness: ", first(ice_model.ice_thickness))
end

set_theme!(Theme(fontsize=24, linewidth=4))
fig = Figure()
axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Ocean temperature (ᵒC)")
axh = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")

lines!(axT, t ./ day, To)
lines!(axh, t ./ day, hi)

display(fig)


#=
# Accumulate data
timeseries = []

function accumulate_timeseries(sim)
    T = model.top_temperature
    h = model.ice_thickness
    push!(timeseries, (time(sim), first(h), first(T)))
end

simulation.callbacks[:save] = Callback(accumulate_timeseries)

run!(simulation)

# Extract and visualize data

t = [datum[1] for datum in timeseries]
h = [datum[2] for datum in timeseries]
T = [datum[3] for datum in timeseries]

set_theme!(Theme(fontsize=24, linewidth=4))

fig = Figure(resolution=(1000, 800))

axT = Axis(fig[1, 1], xlabel="Time (days)", ylabel="Top temperature (ᵒC)")
axh = Axis(fig[2, 1], xlabel="Time (days)", ylabel="Ice thickness (m)")

lines!(axT, t / day, T)
lines!(axh, t / day, h)

display(fig)
=#
