using Oceananigans
using OrthogonalSphericalShellGrids
using ClimaOcean
using ClimaOcean.ECCO
using ClimaOcean.JRA55
using ClimaSeaIce
using ClimaSeaIce.Rheologies
using ClimaSeaIce.SeaIceMomentumEquations   
using Oceananigans.Units
using Printf

using CUDA
CUDA.device!(2)

include("prescribed_external_stress.jl")

arch = GPU()

sea_ice_grid  = TripolarGrid(arch; size=(680, 220, 1), southernmost_latitude=55, z=(-30, 0), halo=(5, 5, 4))
bottom_height = regrid_bathymetry(sea_ice_grid, interpolation_passes=1, minimum_depth=0, major_basins=Inf)
grid = ImmersedBoundaryGrid(sea_ice_grid, GridFittedBottom(bottom_height))

#####
##### Prescribed ocean velocities from ECCO 
#####

# dates = ECCO.all_ECCO_dates(ECCO2Daily())[1:10] # 10 days of evolution

# u_velocity = ECCO.ECCOMetadata(:u_velocity; version=ECCO2Daily(), dates)
# v_velocity = ECCO.ECCOMetadata(:u_velocity; version=ECCO2Daily(), dates)

# uₒ = ECCO.ECCOFieldTimeSeries(u_velocity, arch)
# vₒ = ECCO.ECCOFieldTimeSeries(v_velocity, arch)

ρₒ = 1025.0 # kg/m³
Cᴰ = 5.5e-3 # N/m²

τᵤₒ = τᵥₒ = PrescribedOceanStress(ρₒ, Cᴰ)

##### 
##### Prescribed JRA55 Atmosphere
#####

uₐ = JRA55.JRA55_field_time_series(:eastward_velocity;  architecture=arch, backend=JRA55NetCDFBackend(10))
vₐ = JRA55.JRA55_field_time_series(:northward_velocity; architecture=arch, backend=JRA55NetCDFBackend(10))

ρₐ = 1.125  # kg/m³
Cᴰ = 1.2e-3 # N/m²

τᵤₐ = τᵥₐ = PrescribedAtmosphereStress(uₐ, vₐ, ρₐ, Cᴰ, uₐ.grid)

#####
##### Set up the model
#####

u_bcs = FieldBoundaryConditions(top=nothing, bottom=nothing)
v_bcs = FieldBoundaryConditions(top=nothing, bottom=nothing)

# We use an elasto-visco-plastic rheology and WENO seventh order 
# for advection of h and ℵ
momentum_equations = SeaIceMomentumEquation(grid; 
                                            coriolis = HydrostaticSphericalCoriolis(),
                                            rheology = ElastoViscoPlasticRheology(),
                                            solver   = SplitExplicitSolver(substeps=120))

# Define the model!
model = SeaIceModel(grid; 
                    top_momentum_stress = (u = τᵤₐ, v = τᵥₐ),
                    bottom_momentum_stress = (u = τᵤₒ, v = τᵥₒ),
                    ice_dynamics = momentum_equations,
                    ice_thermodynamics = nothing, # No thermodynamics here
                    advection = WENO(; order = 7),
                    boundary_conditions = (u = u_bcs, v = v_bcs))

# We start with thickenss and concentration from climatology
set!(model.ice_thickness,     ECCOMetadata(:sea_ice_thickness),     inpainting=ClimaOcean.DataWrangling.NearestNeighborInpainting(1))
set!(model.ice_concentration, ECCOMetadata(:sea_ice_area_fraction), inpainting=ClimaOcean.DataWrangling.NearestNeighborInpainting(1))

simulation = Simulation(model, Δt=120, stop_iteration=1) #stop_time=2days)

# # Container to hold the data
htimeseries = []
ℵtimeseries = []
utimeseries = []
vtimeseries = []

## Callback function to collect the data from the `sim`ulation
function accumulate_timeseries(sim)
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    u = sim.model.velocities.u
    v = sim.model.velocities.v
    push!(htimeseries, deepcopy(Array(interior(h))))
    push!(ℵtimeseries, deepcopy(Array(interior(ℵ))))
    push!(utimeseries, deepcopy(Array(interior(u))))
    push!(vtimeseries, deepcopy(Array(interior(v))))
end

wall_time = [time_ns()]

function progress(sim) 
    h = sim.model.ice_thickness
    ℵ = sim.model.ice_concentration
    u = sim.model.velocities.u
    v = sim.model.velocities.v

    hmax = maximum(interior(h))
    ℵmin = minimum(interior(ℵ))
    umax = maximum(interior(u)), maximum(interior(v))
    step_time = 1e-9 * (time_ns() - wall_time[1])

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e), max(h): %.2f, min(ℵ): %.2f, wtime: %s \n",
                   prettytime(sim.model.clock.time),
                   sim.model.clock.iteration,
                   prettytime(sim.Δt),
                   umax..., hmax, ℵmin, prettytime(step_time))

     wall_time[1] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(5))
simulation.callbacks[:save]     = Callback(accumulate_timeseries, IterationInterval(5))

run!(simulation)
