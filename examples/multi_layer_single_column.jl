using Oceananigans
using Oceananigans.Units
using ClimaSeaIce: PrognosticThicknessThermodynamicSeaIceModel, time_step!
using GLMakie

# Create a single column grid with no grid points in the vertical,
# appropriate for a zero heat capacity model.
Nx = Ny = 1
Nz = 20
grid = RectilinearGrid(size = (Nx, Ny, Nz),
                       x = (0, 1),
                       y = (0, 1),
                       z = (0, 1),
                       topology=(Periodic, Periodic, Bounded))

ice_atmosphere_flux = Field{Center, Center, Nothing}(grid)
enthalpy_top_bc = FluxBoundaryCondition(ice_atmosphere_flux)
enthalpy_bcs = FieldBoundaryConditions(top=enthalpy_top_bc)

#bottom_ice_temperature = Field{Center, Center, Nothing}(grid)
temperature_bottom_bc = ValueBoundaryCondition(0)
temperature_bcs = FieldBoundaryConditions(bottom=temperature_bottom_bc)

ρi = 917.0
ci = 2112.2

model = PrognosticThicknessThermodynamicSeaIceModel(grid,
                                                    ice_density = ρi,
                                                    ice_heat_capacity = ci,
                                                    boundary_conditions = (e=enthalpy_bcs, T=temperature_bcs),
                                                    ice_conductivity = 2.0)

z₀ = 1/2
δz = 1/8
Tᵢ(x, y, z) = - 10 * z
set!(model, T=Tᵢ)

Tᵢ(x, y, z) = - 10 * z
set!(model, T=Tᵢ, h=1)

simulation = Simulation(model, Δt=10minute, stop_iteration=100)

T = model.ice_temperature
e = model.ice_enthalpy
h = model.ice_thickness
z = znodes(T)

Nz = size(T, 3)

function compute_ice_atmosphere_flux!(sim)
    T = sim.model.ice_temperature
    Ts = T[1, 1, Nz] + 273.15 # Kelvin
    σ = 5.67e-8 # Stefan-Boltzmann constant
    Q_emission = σ * Ts^4
    Q_total = -500 + Q_emission
    set!(ice_atmosphere_flux, Q_total)
    return nothing
end

# Integrate forward for 1 hour
Tt = []
et = []
ht = Float64[]
tt = Float64[]

function collect_state(sim)
    T = sim.model.ice_temperature
    e = sim.model.ice_enthalpy
    h = sim.model.ice_thickness

    push!(Tt, deepcopy(interior(T, 1, 1, :)))
    push!(et, deepcopy(interior(e, 1, 1, :)))
    push!(ht, deepcopy(first(interior(h))))
    push!(tt, time(sim))

    return nothing
end

simulation.callbacks[:atmos] = Callback(compute_ice_atmosphere_flux!)
simulation.callbacks[:collect] = Callback(collect_state)

progress(sim) = @info string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

run!(simulation)

fig = Figure()
axT = Axis(fig[1, 1], title="T")
axe = Axis(fig[1, 2], title="e")
axh = Axis(fig[2, 1:2], title="h")

Nt = length(tt)

slider = Slider(fig[3, 1:2], range=1:Nt, startvalue=1)
n = slider.value

Tn = @lift Tt[$n]
en = @lift et[$n]
tn = @lift tt[$n]
hn = @lift ht[$n]

lines!(axT, Tn, z)
lines!(axe, en, z)

lines!(axh, tt, ht)
scatter!(axh, tn, hn, marker=:circle, color=:dodgerblue, size=20)

display(fig)

