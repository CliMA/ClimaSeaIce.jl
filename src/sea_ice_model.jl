using Oceananigans.Fields: TracerFields
using Oceananigans.TimeSteppers: TimeStepper
using Oceananigans: tupleit, tracernames
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions
using Oceananigans.Fields: ConstantField
using ClimaSeaIce.SeaIceThermodynamics: external_top_heat_flux
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: flux_summary

struct SeaIceModel{GR, TD, D, TS, CL, U, T, IT, IC, ID, UO, GA, EO, DO, CO, STF, SMS, A} <: AbstractModel{TS}
    grid :: GR
    clock :: CL
    # Prognostic State
    velocities :: U
    tracers :: T
    ice_thickness :: IT
    ice_concentration :: IC
    ice_density :: ID
    # Thermodynamics
    ice_thermodynamics :: TD
    # Dynamics
    ice_dynamics :: D
    ocean_velocities :: UO
    gravitational_acceleration :: GA
    ocean_free_surface :: EO
    ocean_density :: DO
    coriolis :: CO
    # External boundary conditions
    external_heat_fluxes :: STF
    external_momentum_stresses :: SMS
    # Numerics
    timestepper :: TS
    advection :: A
end

function SeaIceModel(grid;
                     clock                      = Clock{eltype(grid)}(time = 0),
                     ice_thickness              = Field{Center, Center, Nothing}(grid),
                     ice_concentration          = Field{Center, Center, Nothing}(grid),
                     ice_salinity               = 0, # psu
                     ice_density                = 900, # kg/m³
                     top_heat_flux              = nothing,
                     bottom_heat_flux           = 0,
                     timestepper                = :RungeKutta3,
                     velocities                 = nothing,
                     advection                  = nothing,
                     top_u_stress               = Field{Face, Center, Nothing}(grid),
                     top_v_stress               = Field{Center, Face, Nothing}(grid),
                     ocean_velocities           = (u = ZeroField(eltype(grid)), v = ZeroField(eltype(grid))),
                     ocean_free_surface         = ZeroField(eltype(grid)),
                     ocean_density              = 1025, # kg/m³
                     coriolis                   = nothing,
                     tracers                    = (),
                     boundary_conditions        = NamedTuple(),
                     gravitational_acceleration = Oceananigans.BuoyancyFormulations.g_Earth,
                     ice_thermodynamics         = SlabSeaIceThermodynamics(grid),
                     ice_dynamics               = SplitExplicitDynamics(grid))

    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)

    # Next, we form a list of default boundary conditions:
    prognostic_field_names = (:u, :v, tracernames(tracers)...)
    default_boundary_conditions = NamedTuple{prognostic_field_names}(Tuple(FieldBoundaryConditions()
                                                                     for name in prognostic_field_names))

    # Then we merge specified, embedded, and default boundary conditions. Specified boundary conditions
    # have precedence, followed by embedded, followed by default.
    boundary_conditions = merge(default_boundary_conditions, boundary_conditions)
    boundary_conditions = regularize_field_boundary_conditions(boundary_conditions, grid, prognostic_field_names)

    if isnothing(velocities) 
        u = Field{Face, Center, Nothing}(grid, boundary_conditions=boundary_conditions.u)
        v = Field{Center, Face, Nothing}(grid, boundary_conditions=boundary_conditions.v)
        velocities = (; u, v)
    end

    tracers = TracerFields(tracers, grid, boundary_conditions)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_salinity in a field 
    ice_salinity = field((Center, Center, Nothing), ice_salinity, grid)

    # Adding thickness and concentration if not there
    prognostic_tracers = merge(tracers, (; h = ice_thickness, ℵ = ice_concentration))
    prognostic_tracers = if ice_salinity isa ConstantField 
        prognostic_tracers 
    else
        merge(prognostic_tracers, (; S = ice_salinity))
    end
    
    # TODO: should we have ice thickness and concentration as part of the tracers or
    # just additional fields of the sea ice model?
    tracers = merge(tracers, (; S = ice_salinity))
    timestepper = TimeStepper(timestepper, grid, prognostic_tracers)

    top_heat_flux = external_top_heat_flux(ice_thermodynamics, top_heat_flux)

    # Package the external fluxes and boundary conditions
    external_heat_fluxes = (top = top_heat_flux,    
                            bottom = bottom_heat_flux) 

    external_momentum_stress = (u = top_u_stress,
                                v = top_v_stress)

    return SeaIceModel(grid,
                       clock,
                       velocities,
                       tracers,
                       ice_thickness,
                       ice_concentration,
                       ice_density,
                       ice_thermodynamics,
                       ice_dynamics,
                       ocean_velocities,
                       convert(eltype(grid), gravitational_acceleration),
                       ocean_free_surface,
                       ocean_density,
                       coriolis,
                       external_heat_fluxes,
                       external_momentum_stress,
                       timestepper,
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
    print(io, "├── ice thermodynamics: ", summary(model.ice_thermodynamics), '\n')
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
