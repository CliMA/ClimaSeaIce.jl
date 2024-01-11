using Oceananigans
using Oceananigans.Architectures: arch_array
using Oceananigans.Fields: ZeroField, ConstantField, VelocityFields
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Units
using Oceananigans.Utils: prettysummary

using SeawaterPolynomials: TEOS10EquationOfState, haline_contraction

using ClimaSeaIce
using ClimaSeaIce: melting_temperature
using ClimaSeaIce.HeatBoundaryConditions: RadiativeEmission, IceWaterThermalEquilibrium
using ClimaSeaIce.SlabSeaIceModels.SlabSeaIceDynamics: FreeDriftRheology

using Printf
using GLMakie
using Statistics

import Oceananigans.ImmersedBoundaries: mask_immersed_field!

mask_immersed_field!(::ConstantField) = nothing

arch = CPU()
Nx = Ny = 256
Lz = 400
x = y = (-50kilometers, 50kilometers)
halo = (4, 4, 4)
topology = (Periodic, Bounded, Bounded)

ice_grid = RectilinearGrid(arch; x, y,
                           size = (Nx, Ny),
                           topology = (topology[1], topology[2], Flat),
                           halo = halo[1:2])

mask(x, y, z) = abs(x) < 5kilometers && y < -10kilometers

ice_grid = ImmersedBoundaryGrid(ice_grid, GridFittedBoundary(mask))

coriolis = FPlane(f=1.4e-4)

u, v, w = VelocityFields(ice_grid)
ocean_surface_velocities = (; u, v)

# Ocean pushing at 0.1 m/s in the zonal direction (westward)
fill!(u, 0.1)

ice_model = SlabSeaIceModel(ice_grid;
                            ocean_velocities = ocean_surface_velocities,
                            advection = WENO(),
                            rheology = FreeDriftRheology(10), # Ice is just drifting
                            internal_heat_flux = ConductiveFlux(conductivity=2),
                            top_heat_flux = ConstantField(0), # W m⁻²
                            top_heat_boundary_condition = PrescribedTemperature(0),
                            bottom_heat_boundary_condition = PrescribedTemperature(0),
                            bottom_heat_flux = ConstantField(0),
                            coriolis = nothing)

ice_simulation = Simulation(ice_model, Δt=20minutes, stop_time = 1days)

hᵢ(x, y) = 2.0
ℵᵢ(x, y) = 1.0

set!(ice_model, h=hᵢ, ℵ=ℵᵢ)

h = ice_model.ice_thickness
ℵ = ice_model.concentration

function progress(sim)
    h = sim.model.ice_thickness
    ℵ = sim.model.concentration
    u = sim.model.velocities.u
    v = sim.model.velocities.v
    msg1 = @sprintf("Iter: % 6d, time: % 12s", iteration(sim), prettytime(sim))
    msg2 = @sprintf(", max(h): %.2f", maximum(h))
    msg3 = @sprintf(", max(ℵ): %.2f", maximum(ℵ))
    msg6 = @sprintf(", maximum(vel): (%.2f, %.2f)", maximum(abs, u), maximum(abs, v))
    @info msg1 * msg2 * msg3 * msg6 
    return nothing
end

ice_simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))

h = ice_model.ice_thickness
ℵ = ice_model.concentration
ui, vi = ice_model.velocities

ht = []
ℵt = []
ut = []
vt = []
tt = []

function saveoutput(sim)
    hn = Array(interior(h, :, :, 1))
    ℵn = Array(interior(ℵ, :, :, 1))
    unᵢ = Array(interior(ui, :, :, 1))
    vnᵢ = Array(interior(vi, :, :, 1))
    push!(ht, hn)
    push!(ℵt, ℵn)
    push!(ut, unᵢ)
    push!(vt, vnᵢ)
    push!(tt, time(sim))
end

ice_simulation.callbacks[:output] = Callback(saveoutput, IterationInterval(1))

run!(ice_simulation)

#####
##### Viz
#####

set_theme!(Theme(fontsize=24))

x = xnodes(ice_grid, Center())
y = ynodes(ice_grid, Center())

fig = Figure(resolution=(2400, 700))

axh = Axis(fig[1, 1], xlabel="x (km)", ylabel="y (km)", title="Ice thickness")
axℵ = Axis(fig[1, 2], xlabel="x (km)", ylabel="y (km)", title="Ice Concentration")
axU = Axis(fig[1, 3], xlabel="x (km)", ylabel="y (km)", title="Ice zonal velocity")
axV = Axis(fig[1, 4], xlabel="x (km)", ylabel="y (km)", title="Ice meridional velocity")

Nt = length(tt)
slider = Slider(fig[2, 1:4], range=1:Nt, startvalue=Nt)
n = slider.value

title = @lift string("Melt-driven baroclinic instability after ", prettytime(tt[$n]))
Label(fig[0, 1:3], title)

hn = @lift ht[$n][:, 50]
ℵn = @lift ℵt[$n][:, 50]
un = @lift ut[$n][:, 50]
vn = @lift vt[$n]
Un = @lift mean(ut[$n], dims=1)[:]

x = x ./ 1e3
y = y ./ 1e3

# heatmap!(axh, x, y, hn, colorrange=(0, 5), colormap=:grays)
lines!(axh, x, hn) #, colorrange=(0, 1), colormap=:grays)
ylims!(axh, (0, 5))
lines!(axℵ, x, ℵn) #, colorrange=(0, 1), colormap=:grays)
lines!(axU, x, un) #, colorrange=(0, 1), colormap=:grays)
# heatmap!(axℵ, x, y, ℵn, colorrange=(0, 1), colormap=:grays)
# heatmap!(axU, x, y, un, colorrange = (-0.1, 0.1), colormap=:haline)
heatmap!(axV, x, y, vn, colorrange = (-0.1, 0.1), colormap=:haline)

display(fig)

#=
record(fig, "salty_baroclinic_ice_cube.mp4", 1:Nt, framerate=48) do nn
    @info string(nn)
    n[] = nn
end
=#

