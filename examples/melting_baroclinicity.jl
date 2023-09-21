using Oceananigans
using Oceananigans.Units
using Oceananigans.Utils: prettysummary
using Oceananigans.Fields: ZeroField
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Printf
using ClimaSeaIce
using ClimaSeaIce: melting_temperature
using SeawaterPolynomials: TEOS10EquationOfState
using ClimaSeaIce.ThermalBoundaryConditions: RadiativeEmission, IceWaterThermalEquilibrium
using GLMakie

import Oceananigans.Simulations: time_step!, time

include("ice_ocean_model.jl")

Nx = 64
Ny = 64
x = (0, 100kilometers)
y = (-50kilometers, 50kilometers)
Δt = 4hours
mixed_layer_depth = hₒ = 5
halo = (4, 4, 4)
topology = (Periodic, Bounded, Bounded)
ice_topology = (topology[1], topology[2], Flat)

ice_grid = RectilinearGrid(size=(Nx, Ny); x, y, topology = ice_topology, halo = halo[1:2])
ocean_grid = RectilinearGrid(size=(Nx, Ny, 1); x, y, z=(-hₒ, 0), topology, halo)

# Top boundary conditions:
#   - outgoing radiative fluxes emitted from surface
#   - incoming shortwave radiation starting after 40 days

radiative_emission = RadiativeEmission()
ice_ocean_heat_flux = Field{Center, Center, Nothing}(ice_grid)
top_ocean_heat_flux = Qᵀ = Field{Center, Center, Nothing}(ice_grid)
top_salt_flux = Qˢ = Field{Center, Center, Nothing}(ice_grid)
solar_insolation = I₀ = Field{Center, Center, Nothing}(ice_grid)
set!(solar_insolation, -900) # W m⁻²

# Generate a zero-dimensional grid for a single column slab model 

boundary_conditions = (T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                       u = FieldBoundaryConditions(top=FluxBoundaryCondition(-1e-6)),
                       S = FieldBoundaryConditions(top=FluxBoundaryCondition(Qˢ)))

equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

ocean_model = HydrostaticFreeSurfaceModel(grid = ocean_grid; buoyancy, boundary_conditions,
                                          momentum_advection = WENO(),
                                          tracer_advection = WENO(),
                                          closure = nothing,
                                          coriolis = FPlane(f=-1e-4),
                                          tracers = (:T, :S, :e))

Nz = size(ocean_grid, 3)
So = ocean_model.tracers.S
ocean_surface_salinity = view(So, :, :, Nz)
bottom_bc = IceWaterThermalEquilibrium(ocean_surface_salinity)

u, v, w = ocean_model.velocities
ocean_surface_velocities = (u = view(u, :, :, Nz),
                            v = view(v, :, :, Nz),    
                            w = ZeroField())

ice_model = SlabSeaIceModel(ice_grid;
                            velocities = ocean_surface_velocities,
                            advection = WENO(),
                            ice_consolidation_thickness = 0.1,
                            ice_salinity = 0,
                            internal_thermal_flux = ConductiveFlux(conductivity=2),
                            top_thermal_flux = (solar_insolation, radiative_emission),
                            bottom_thermal_boundary_condition = bottom_bc,
                            bottom_thermal_flux = ice_ocean_heat_flux)

ocean_simulation = Simulation(ocean_model; Δt=20minutes, verbose=false)
ice_simulation = Simulation(ice_model, Δt=20minutes, verbose=false)


# Initial condition
S₀ = 30
T₀ = melting_temperature(ice_model.phase_transitions.liquidus, S₀)

uᵢ(x, y, z) = 0.1 * randn()
vᵢ(x, y, z) = 0.1 * randn()
Tᵢ(x, y, z) = T₀ + 0.1 * randn()
Sᵢ(x, y, z) = S₀ + 0.1 * randn()
hᵢ(x, y) = y < 0 ? 3 : 0

set!(ocean_model, S=Sᵢ, T=T₀)
set!(ice_model, h=hᵢ)

coupled_model = IceOceanModel(ice_simulation, ocean_simulation)
coupled_simulation = Simulation(coupled_model, Δt=10minutes, stop_time=10day)

using SeawaterPolynomials: thermal_expansion, haline_contraction

S = ocean_model.tracers.S
β = haline_contraction(T₀, S₀, 0, equation_of_state)
g = ocean_model.buoyancy.model.gravitational_acceleration
by = - g * β * ∂y(S)

function progress(sim)
    h = sim.model.ice.model.ice_thickness
    S = sim.model.ocean.model.tracers.S
    msg1 = @sprintf("Iter: % 6d, time: % 12s", iteration(sim), prettytime(sim))
    msg2 = @sprintf(", max(h): %.2f", maximum(h))
    msg3 = @sprintf(", min(S): %.2f", minimum(S))
    msg4 = @sprintf(", max|∂y b|: %.2e", maximum(abs, by))
    @info msg1 * msg2 * msg3 * msg4
    return nothing
end

coupled_simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

h = ice_model.ice_thickness
T = ocean_model.tracers.T
S = ocean_model.tracers.S
u, v, w = ocean_model.velocities
η = ocean_model.free_surface.η

ht = []
Tt = []
St = []
ut = []
vt = []
ηt = []
tt = []

function saveoutput(sim)
    hn = deepcopy(interior(h, :, :, 1))
    Tn = deepcopy(interior(T, :, :, 1))
    Sn = deepcopy(interior(S, :, :, 1))
    un = deepcopy(interior(u, :, :, 1))
    vn = deepcopy(interior(v, :, :, 1))
    ηn = deepcopy(interior(η, :, :, 1))
    push!(ht, hn)
    push!(Tt, Tn)
    push!(St, Sn)
    push!(ut, un)
    push!(vt, vn)
    push!(ηt, ηn)
    push!(tt, time(sim))
end

coupled_simulation.callbacks[:output] = Callback(saveoutput, IterationInterval(10))

run!(coupled_simulation)

using GLMakie

fig = Figure(resolution=(1800, 900))

axh = Axis(fig[1, 1])
axT = Axis(fig[1, 2])
axS = Axis(fig[1, 3])
axu = Axis(fig[2, 1])
axv = Axis(fig[2, 2])
axη = Axis(fig[2, 3])


Nt = length(tt)
slider = Slider(fig[3, 1:3], range=1:Nt, startvalue=1)
n = slider.value

hn = @lift ht[$n]
Tn = @lift Tt[$n]
Sn = @lift St[$n]
un = @lift ut[$n]
vn = @lift vt[$n]
ηn = @lift ηt[$n]

heatmap!(axh, hn)
heatmap!(axT, Tn)
heatmap!(axS, Sn)
heatmap!(axu, un)
heatmap!(axv, vn)
heatmap!(axη, ηn)

display(fig)

