#=
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

Nx = Ny = 128
x = (0, 50kilometers)
y = (-25kilometers, 25kilometers)
Δt = 4hours
halo = (4, 4, 4)
topology = (Periodic, Bounded, Bounded)
ice_topology = (topology[1], topology[2], Flat)

ice_grid = RectilinearGrid(size=(Nx, Ny); x, y, topology = ice_topology, halo = halo[1:2])
ocean_grid = RectilinearGrid(size=(Nx, Ny, 2);  topology, halo, x, y,
                             z=[-6, -2, 0])

# Top boundary conditions:
#   - outgoing radiative fluxes emitted from surface
#   - incoming shortwave radiation starting after 40 days

radiative_emission = RadiativeEmission()
ice_ocean_heat_flux      = Field{Center, Center, Nothing}(ice_grid)
top_ocean_heat_flux = Qᵀ = Field{Center, Center, Nothing}(ice_grid)
top_salt_flux       = Qˢ = Field{Center, Center, Nothing}(ice_grid)
solar_insolation    = I₀ = Field{Center, Center, Nothing}(ice_grid)
set!(solar_insolation, -600) # W m⁻²

# Generate a zero-dimensional grid for a single column slab model 

boundary_conditions = (T = FieldBoundaryConditions(top=FluxBoundaryCondition(Qᵀ)),
                       S = FieldBoundaryConditions(top=FluxBoundaryCondition(Qˢ)))

equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

ocean_model = HydrostaticFreeSurfaceModel(grid = ocean_grid; buoyancy, boundary_conditions,
                                          momentum_advection = WENO(),
                                          tracer_advection = WENO(),
                                          closure = nothing,
                                          coriolis = FPlane(f=1.4e-4),
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
                            ice_consolidation_thickness = 0.05,
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

Tᵢ(x, y, z) = T₀ + 0.1 * randn()
Sᵢ(x, y, z) = S₀ + 0.1 * randn()
hᵢ(x, y) = y < 0 ? 3 : 0

set!(ocean_model, S=Sᵢ, T=T₀)
set!(ice_model, h=hᵢ)

coupled_model = IceOceanModel(ice_simulation, ocean_simulation)
coupled_simulation = Simulation(coupled_model, Δt=1minutes, stop_time=30days)

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
ζt = []
tt = []

ζ = Field(∂x(v) - ∂y(u))

function saveoutput(sim)
    compute!(ζ)
    hn = deepcopy(interior(h, :, :, 1))
    Tn = deepcopy(interior(T, :, :, 1))
    Sn = deepcopy(interior(S, :, :, 1))
    un = deepcopy(interior(u, :, :, 1))
    vn = deepcopy(interior(v, :, :, 1))
    ηn = deepcopy(interior(η, :, :, 1))
    ζn = deepcopy(interior(ζ, :, :, 1))
    push!(ht, hn)
    push!(Tt, Tn)
    push!(St, Sn)
    push!(ut, un)
    push!(vt, vn)
    push!(ηt, ηn)
    push!(ζt, ζn)
    push!(tt, time(sim))
end

coupled_simulation.callbacks[:output] = Callback(saveoutput, IterationInterval(10))

run!(coupled_simulation)
=#

using GLMakie
using Statistics

set_theme!(Theme(fontsize=24))

x = xnodes(ocean_grid, Center())
y = ynodes(ocean_grid, Center())

fig = Figure(resolution=(2400, 800))

axh = Axis(fig[1, 1], xlabel="x (km)", ylabel="y (km)", title="Ice thickness")
axT = Axis(fig[1, 2], xlabel="x (km)", ylabel="y (km)", title="Ocean temperature")
axS = Axis(fig[1, 3], xlabel="x (km)", ylabel="y (km)", title="Ocean salinity")
#axZ = Axis(fig[1, 4])

#axU = Axis(fig[2, 0])
# axu = Axis(fig[2, 1])
# axv = Axis(fig[2, 2])
# axη = Axis(fig[2, 3])

Nt = length(tt)
slider = Slider(fig[2, 1:3], range=1:Nt, startvalue=1)
n = slider.value

title = @lift string("Melt-driven baroclinic instability after ", prettytime(tt[$n]))
Label(fig[0, 1:3], title)

hn = @lift ht[$n]
Tn = @lift Tt[$n]
Sn = @lift St[$n]
un = @lift ut[$n]
vn = @lift vt[$n]
ηn = @lift ηt[$n]
ζn = @lift ζt[$n]
Un = @lift mean(ut[$n], dims=1)[:]

x = x ./ 1e3
y = y ./ 1e3

heatmap!(axh, x, y, hn, colorrange=(0, 1), colormap=:grays)
heatmap!(axT, x, y, Tn, colorrange=(-1.7, -1.5), colormap=:thermal)
heatmap!(axS, x, y, Sn, colorrange=(28, 30), colormap=:haline)
#heatmap!(axZ, ζn, colormap=:redblue)

#lines!(axU, Un, y)
#heatmap!(axu, un)
#heatmap!(axv, vn)
#heatmap!(axη, ηn)

display(fig)

record(fig, "melting_baroclinicity.mp4", 1:Nt, framerate=48) do nn
    @info string(nn)
    n[] = nn
end

