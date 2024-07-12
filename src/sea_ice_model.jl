using Oceananigans.Fields: TracerFields
using ClimaSeaIce.SeaIceThermodynamics: external_top_heat_flux
using Oceananigans: tupleit, tracernames
using Oceananigans.TimeSteppers: TimeStepper

struct SeaIceModel{GR, TD, D, CL, TS, U, T, IT, IC, ID, UO, DO, CO, STF, SMS, A} <: AbstractModel{TS}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # Prognostic State
    velocities :: U
    tracers :: T
    ice_thickness :: IT
    ice_concentration :: IC
    ice_density :: ID
    # Thermodynamics
    sea_ice_thermodynamics :: TD
    # Dynamics
    sea_ice_dynamics :: D
    ocean_velocities :: UO
    ocean_density :: DO
    coriolis :: CO
    # External boundary conditions
    external_heat_fluxes :: STF
    external_momentum_stresses :: SMS
    # Numerics
    advection :: A
end

function SeaIceModel(grid;
                     clock                  = Clock{eltype(grid)}(time = 0),
                     ice_thickness          = Field{Center, Center, Nothing}(grid),
                     ice_concentration      = Field{Center, Center, Nothing}(grid),
                     ice_salinity           = 0, # psu
                     ice_density            = 900, # kg/m³
                     top_heat_flux          = nothing,
                     bottom_heat_flux       = 0,
                     velocities             = nothing,
                     advection              = nothing,
                     top_u_stress           = Field{Face, Center, Nothing}(grid),
                     top_v_stress           = Field{Center, Face, Nothing}(grid),
                     ocean_velocities       = (u = ZeroField(), v = ZeroField()),
                     ocean_density          = 1025, # kg/m³
                     coriolis               = nothing,
                     tracers                = (),
                     boundary_conditions    = NamedTuple(),
                     sea_ice_thermodynamics = SlabSeaIceThermodynamics(grid),
                     sea_ice_dynamics       = ExplicitMomentumSolver(grid))

    if isnothing(velocities) 
        u = XFaceField(grid)
        v = YFaceField(grid)
        velocities = (; u, v)
    end

    # Only one time-stepper is supported currently
    timestepper = ForwardEulerTimestepper()

    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)
    tracers = TracerFields(tracers, grid, boundary_conditions)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_salinity in a field
    ice_salinity = field((Center, Center, Nothing), ice_salinity, grid)

    # Adding thickness and concentration if not there
    tracer_names = tuple(unique((tracernames(tracers)..., :h, :ℵ))...)

    # Only one time-stepper is supported currently
    timestepper = TimeStepper(:QuasiAdamsBashforth2, grid, tracer_names;
                              Gⁿ = TracerFields(tracer_names, grid),
                              G⁻ = TracerFields(tracer_names, grid))

    top_heat_flux = external_top_heat_flux(sea_ice_thermodynamics, top_heat_flux)

    # Package the external fluxes and boundary conditions
    external_heat_fluxes = (top = top_heat_flux,    
                            bottom = bottom_heat_flux) 

    external_momentum_stress = (u = top_u_stress,
                                v = top_v_stress)

    return SeaIceModel(grid,
                       clock,
                       timestepper,
                       velocities,
                       tracers,
                       ice_thickness,
                       ice_concentration,
                       ice_density,
                       sea_ice_thermodynamics,
                       sea_ice_dynamics,
                       ocean_velocities,
                       ocean_density,
                       coriolis,
                       external_heat_fluxes,
                       external_momentum_stress,
                       advection)
end

const SIM = SeaIceModel

function set!(model::SIM; h=nothing, ℵ=nothing)
    !isnothing(h) && set!(model.ice_thickness, h)
    !isnothing(ℵ) && set!(model.ice_concentration, ℵ)
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
    print(io, "├── sea ice thermodynamics: ", summary(model.sea_ice_thermodynamics), '\n')
    print(io, "├── sea ice dynamics: ", summary(model.sea_ice_dynamics), '\n')
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
                           fields(model.sea_ice_thermodynamics))

# TODO: make this correct
prognostic_fields(model::SIM) = merge((; h  = model.ice_thickness,
                                         ℵ  = model.ice_concentration),
                                      model.tracers,
                                      model.velocities)
