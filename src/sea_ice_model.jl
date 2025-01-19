using Oceananigans.Fields: TracerFields
using Oceananigans.TimeSteppers: TimeStepper
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions
using ClimaSeaIce.SeaIceThermodynamics: PrescribedTemperature
using Oceananigans: tupleit, tracernames
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions
using Oceananigans.Forcings: model_forcing
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: flux_summary

struct SeaIceModel{GR, TD, D, TS, CL, U, T, IT, IC, ID, STF, SMS, A, F} <: AbstractModel{TS}
    grid :: GR
    clock :: CL
    forcing :: F
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
    # External boundary conditions
    external_heat_fluxes :: STF
    external_momentum_stresses :: SMS
    # Numerics
    timestepper :: TS
    advection :: A
end

function SeaIceModel(grid;
                     clock                  = Clock{eltype(grid)}(time = 0),
                     ice_thickness          = nothing,
                     ice_concentration      = nothing,
                     ice_salinity           = 0, # psu
                     ice_density            = 900, # kg m⁻³
                     top_heat_flux          = nothing,
                     bottom_heat_flux       = 0,
                     velocities             = nothing,
                     timestepper            = :QuasiAdamsBashforth2,
                     advection              = nothing,
                     top_momentum_stress    = (u = nothing, v = nothing),
                     bottom_momentum_stress = (u = nothing, v = nothing),
                     tracers                = (),
                     boundary_conditions    = NamedTuple(),
                     ice_thermodynamics     = SlabSeaIceThermodynamics(grid),
                     ice_dynamics           = nothing,
                     forcing                = NamedTuple())

    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)

    # Next, we form a list of default boundary conditions:
    field_names = (:u, :v, :h, :ℵ, :S, tracernames(tracers)...)
    default_boundary_conditions = NamedTuple{field_names}(Tuple(FieldBoundaryConditions()
                                                          for name in field_names))

    # Then we merge specified, embedded, and default boundary conditions. Specified boundary conditions
    # have precedence, followed by embedded, followed by default.
    boundary_conditions = merge(default_boundary_conditions, boundary_conditions)
    boundary_conditions = regularize_field_boundary_conditions(boundary_conditions, grid, field_names)
    
    if isnothing(velocities) 
        u = Field{Face, Center, Center}(grid, boundary_conditions=boundary_conditions.u)
        v = Field{Center, Face, Center}(grid, boundary_conditions=boundary_conditions.v)
        velocities = (; u, v)
    end

    tracers = TracerFields(tracers, grid, boundary_conditions)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_salinity in a field 
    ice_salinity = field((Center, Center, Center), ice_salinity, grid)
    ice_density  = field((Center, Center, Center), ice_density, grid)

    # Construct prognostic fields if not provided
    ice_thickness = isnothing(ice_thickness) ?  Field{Center, Center, Center}(grid, boundary_conditions=boundary_conditions.h) : ice_thickness
    ice_concentration = isnothing(ice_concentration) ? Field{Center, Center, Center}(grid, boundary_conditions=boundary_conditions.ℵ) : ice_concentration

    # Adding thickness and concentration if not there
    prognostic_fields = merge(tracers, (; h = ice_thickness, ℵ = ice_concentration))
    prognostic_fields = if ice_salinity isa ConstantField 
        prognostic_fields 
    else
        merge(prognostic_fields, (; S = ice_salinity))
    end
    
    prognostic_fields = isnothing(ice_dynamics) ? prognostic_fields : merge(prognostic_fields, velocities)

    # TODO: should we have ice thickness and concentration as part of the tracers or
    # just additional fields of the sea ice model?
    tracers = merge(tracers, (; S = ice_salinity))
    timestepper = TimeStepper(timestepper, grid, prognostic_fields)

    if !isnothing(ice_thermodynamics)
        if isnothing(top_heat_flux)
            if ice_thermodynamics.heat_boundary_conditions.top isa PrescribedTemperature
                # Default: external top flux is in equilibrium with internal fluxes
                top_heat_flux = ice_thermodynamics.internal_heat_flux
            else
                # Default: no external top surface flux
                top_heat_flux = 0
            end
        end
    end

    forcing = model_forcing(prognostic_fields; forcing...)
    
    # Package the external fluxes and boundary conditions
    external_heat_fluxes = (top = top_heat_flux,    
                            bottom = bottom_heat_flux) 

    external_momentum_stresses = (top = top_momentum_stress,
                                  bottom = bottom_momentum_stress)

    return SeaIceModel(grid,
                       clock,
                       forcing, 
                       velocities,
                       tracers,
                       ice_thickness,
                       ice_concentration,
                       ice_density,
                       ice_thermodynamics,
                       ice_dynamics,
                       external_heat_fluxes,
                       external_momentum_stresses,
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
    print(io, "├── ice_thermodynamics: ", summary(model.ice_thermodynamics), '\n')
    print(io, "├── advection: ", summary(model.advection), '\n')
    print(io, "└── external_heat_fluxes: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_heat_fluxes.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_heat_fluxes.bottom, "     "))
end
         
reset!(::SIM) = nothing
initialize!(::SIM) = nothing
default_included_properties(::SIM) = tuple(:grid)

# Fallback
fields(::Nothing) = NamedTuple()

fields(model::SIM) = merge((; h  = model.ice_thickness,
                              ℵ  = model.ice_concentration),
                           model.tracers,
                           model.velocities,
                           fields(model.ice_thermodynamics),
                           fields(model.ice_dynamics))

# TODO: make this correct
prognostic_fields(model::SIM) = merge((; h  = model.ice_thickness,
                                         ℵ  = model.ice_concentration),
                                      model.tracers,
                                      model.velocities)
