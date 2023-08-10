using ClimaSeaIce: ThermodynamicSeaIceModel, MolecularDiffusivity, temperature_flux
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators: Δzᶜᶜᶠ
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Grids: znode
using GLMakie

#####
##### Set up an Ocean column model
#####

ocean_grid = RectilinearGrid(size=20, z=(-100, 0), topology=(Flat, Flat, Bounded))
equation_of_state = LinearEquationOfState(thermal_expansion=2e-4, haline_contraction=8e-5)

# Set up boundary conditions ready for coupling (on the CPU)
top_T_flux = zeros(1, 1)
top_S_flux = zeros(1, 1)
T_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(top_T_flux))
S_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(top_S_flux))

ocean_model = HydrostaticFreeSurfaceModel(grid = ocean_grid,
                                          closure = CATKEVerticalDiffusivity(),
                                          tracers = (:T, :S, :e),
                                          buoyancy = SeawaterBuoyancy(; equation_of_state),
                                          boundary_conditions = (; T=T_bcs, S=S_bcs))


# Ocean initial condition
α = equation_of_state.thermal_expansion
g = ocean_model.buoyancy.model.gravitational_acceleration
N² = 1e-6
dTdz = - α * g * N²
T₀ = 0.1 # Temperature at surface
Tᵢ(x, y, z) = T₀ + dTdz * z
Sᵢ(x, y, z) = 35
set!(ocean_model, T=Tᵢ, S=Sᵢ)

#####
##### Set up a ThermodynamicSeaIceModel
#####

# Build a grid with 10 cm resolution
ice_grid = RectilinearGrid(size=20, z=(-1, 0), topology=(Flat, Flat, Bounded))

# Set up a simple problem and build the ice model
ice_closure = MolecularDiffusivity(ice_grid, κ_ice=1e-5, κ_water=1e-6)

#####
##### Create temperature boundary conditions
#####

# Atmopshere-ice surface temperature parameters
initial_air_ice_temperature = -5
air_ice_temperature_fluctuation_amplitude = 5
air_ice_tendency = -0.5 / day # ᵒC s⁻¹

parameters = (Tᵢ = initial_air_ice_temperature,
              Tₐ = air_ice_temperature_fluctuation_amplitude,
              dTdt = air_ice_tendency)

@inline air_ice_temperature(x, y, t, p) = p.dTdt * t + p.Tₐ * sin(2π*t / day) + p.Tᵢ
top_T_bc = ValueBoundaryCondition(air_ice_temperature; parameters)
T_bcs = FieldBoundaryConditions(top=top_T_bc)

ice_model = ThermodynamicSeaIceModel(grid = ice_grid,
                                     closure = ice_closure,
                                     boundary_conditions=(; T=T_bcs))


κ = 1e-5
Δz = Δzᶜᶜᶠ(1, 1, 1, ice_grid)
Δt = 0.1 * Δz^2 / κ

ice_simulation = Simulation(ice_model; Δt, stop_iteration=10)

@inline is_ocean(i, j, k, grid, ϕ) = @inbounds ϕ[i, j, k] == 1

function fix_liquid_temperature!(ϕ_ice_model, T_ice_model, grid, T_ocean_model) 
    Nz = size(grid, 3)

    i = j = 1
    @inbounds for k = 1:Nz
        T_ice_model[i, j, k] = ifelse(is_ocean(i, j, k, grid, ϕ_ice_model),
                                      T_ocean_model[i, j],
                                      T_ice_model[i, j, k])
    end

    return nothing
end

function set_liquid_state!(ice_simulation)
    ice_model = ice_simulation.model

    # Extract ocean temperature at top grid point
    ocean_Nz = size(ocean_model.grid, 3)
    ocean_temperature = interior(ocean_model.tracers.T, :, :, ocean_Nz)

    # Eventually this will launch a kernel on the GPU
    ϕ = ice_model.state.ϕ
    T = ice_model.state.T
    fix_liquid_temperature!(ϕ, T, ice_model.grid, ocean_temperature)


    return nothing
end

ice_simulation.callbacks[:ocean_state] = Callback(set_liquid_state!)

#=
ocean_simulation = Simulation(ocean_model; Δt)
=#
