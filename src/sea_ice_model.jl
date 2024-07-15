using Oceananigans.Fields: TracerFields
using ClimaSeaIce.SeaIceThermodynamics: external_top_heat_flux
using Oceananigans: tupleit

struct SeaIceModel{GR, TD, D, CL, TS, U, T, IT, IC, STF, SMS, A} <: AbstractModel{TS}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # Prognostic State
    velocities :: U
    tracers :: T
    ice_thickness :: IT
    ice_concentration :: IC
    # Thermodynamics
    ice_thermodynamics :: TD
    # Dynamics
    ice_dynamics :: D
    # External boundary conditions
    external_heat_fluxes :: STF
    external_momentum_stresses :: SMS
    # Numerics
    advection :: A
end

function SeaIceModel(grid;
                     clock               = Clock{eltype(grid)}(time = 0),
                     ice_thickness       = Field{Center, Center, Nothing}(grid),
                     ice_concentration   = Field{Center, Center, Nothing}(grid),
                     ice_salinity        = 0, # psu
                     top_heat_flux       = nothing,
                     bottom_heat_flux    = 0,
                     velocities          = nothing,
                     advection           = nothing,
                     top_momentum_stress = nothing, # Fix when introducing dynamics
                     tracers             = (),
                     boundary_conditions = NamedTuple(),
                     ice_thermodynamics  = SlabSeaIceThermodynamics(grid),
                     ice_dynamics        = nothing)

    if isnothing(velocities)
        velocities = (u = ZeroField(), v=ZeroField(), w=ZeroField())
    end

    # Only one time-stepper is supported currently
    timestepper = ForwardEulerTimestepper()

    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)
    tracers = TracerFields(tracers, grid, boundary_conditions)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_salinity in a field
    ice_salinity = field((Center, Center, Nothing), ice_salinity, grid)
    tracers = merge(tracers, (; S = ice_salinity))

    top_heat_flux = external_top_heat_flux(ice_thermodynamics, top_heat_flux)

    # Package the external fluxes and boundary conditions
    external_heat_fluxes = (top = top_heat_flux,    
                            bottom = bottom_heat_flux) 

    return SeaIceModel(grid,
                       clock,
                       timestepper,
                       velocities,
                       tracers,
                       ice_thickness,
                       ice_concentration,
                       ice_thermodynamics,
                       ice_dynamics,
                       external_heat_fluxes,
                       top_momentum_stress,
                       advection)
end

const SIM = SeaIceModel

function set!(model::SIM; h=nothing, ℵ=nothing)
    !isnothing(h) && set!(model.ice_thickness, h)
    !isnothing(ℵ) && set!(model.ice_conentration, ℵ)
    return nothing
end

Base.summary(model::SIM) = "SeaIceModel"
prettytime(model::SIM) = prettytime(model.clock.time)
iteration(model::SIM) = model.clock.iteration

function Base.show(io::IO, model::SIM)
    grid = model.grid
    arch = architecture(grid)
    gridname = typeof(grid).name.wrapper
    timestr = string("(time = ", prettytime(model), ", iteration = ", iteration(model), ")")

    print(io, "SeaIceModel{", typeof(arch), ", ", gridname, "}", timestr, '\n')
    print(io, "├── grid: ", summary(model.grid), '\n')
    print(io, "├── ice thermodynamics: ", summary(model.ice_thermodynamics), '\n')
    print(io, "├── ice dynamics: ", summary(model.ice_dynamics), '\n')
    print(io, "├── advection: ", summary(model.advection), '\n')
    print(io, "└── external_heat_fluxes: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_heat_fluxes.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_heat_fluxes.bottom, "     "))
end
         
reset!(::SIM) = nothing
initialize!(::SIM) = nothing
default_included_properties(::SIM) = tuple(:grid)

fields(model::SIM) = merge((; h  = model.ice_thickness,
                              ℵ  = model.ice_concentration),
                           model.tracers,
                           model.velocities,
                           fields(model.ice_thermodynamics))

# TODO: make this correct
prognostic_fields(model::SIM) = merge((; h  = model.ice_thickness,
                                         ℵ  = model.ice_concentration),
                                      model.tracers,
                                      model.velocities)
