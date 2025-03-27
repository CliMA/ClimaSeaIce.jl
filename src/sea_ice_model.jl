using Oceananigans.Architectures: architecture
using Oceananigans.Fields: TracerFields
using Oceananigans.TimeSteppers: TimeStepper
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions
using ClimaSeaIce.SeaIceThermodynamics: PrescribedTemperature
using Oceananigans: tupleit, tracernames
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions
using Oceananigans.Forcings: model_forcing
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: flux_summary
using ClimaSeaIce.Rheologies: rheology_prognostic_tracers

struct SeaIceModel{GR, TD, D, TS, CL, U, T, IT, IC, ID, CT, STF, A, F, Arch} <: AbstractModel{TS, Arch}
    architecture :: Arch
    grid :: GR
    clock :: CL
    forcing :: F
    # Prognostic State
    velocities :: U
    tracers :: T
    ice_thickness :: IT
    ice_concentration :: IC
    ice_density :: ID
    ice_consolidation_thickness :: CT
    # Thermodynamics
    ice_thermodynamics :: TD
    # Dynamics
    dynamics :: D
    # External boundary conditions
    external_heat_fluxes :: STF
    # Numerics
    timestepper :: TS
    advection :: A
end

assumed_sea_ice_field_location(name) = name === :u  ? (Face,   Face,   Nothing) :
                                       name === :v  ? (Face,   Face,   Nothing) :
                                                      (Center, Center, Nothing)

function SeaIceModel(grid;
                     clock                       = Clock{eltype(grid)}(time = 0),
                     ice_consolidation_thickness = 0.05, # m
                     ice_salinity                = 0, # psu
                     ice_density                 = 900, # kg m⁻³
                     top_heat_flux               = nothing,
                     bottom_heat_flux            = 0,
                     velocities                  = nothing,
                     advection                   = nothing,
                     tracers                     = (),
                     boundary_conditions         = NamedTuple(),
                     ice_thermodynamics          = SlabSeaIceThermodynamics(grid),
                     dynamics                    = nothing,
                     forcing                     = NamedTuple())

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_consolidation_thickness in a field
    ice_consolidation_thickness = field((Center, Center, Nothing), ice_consolidation_thickness, grid)

    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)
    tracers = (tracers..., rheology_prognostic_tracers(dynamics)...) # add prognostic tracers

    # Next, we form a list of default boundary conditions:
    field_names = (:u, :v, :h, :ℵ, :S, tracernames(tracers)...)

    bc_tuple = Tuple(FieldBoundaryConditions(grid, assumed_sea_ice_field_location(name))
                     for name in field_names)
    default_boundary_conditions = NamedTuple{field_names}(bc_tuple)

    # Then we merge specified, embedded, and default boundary conditions. Specified boundary conditions
    # have precedence, followed by embedded, followed by default.
    boundary_conditions = merge(default_boundary_conditions, boundary_conditions)
    boundary_conditions = regularize_field_boundary_conditions(boundary_conditions, grid, field_names)
    
    if isnothing(velocities) 
        u = Field{Face, Face, Nothing}(grid, boundary_conditions=boundary_conditions.u)
        v = Field{Face, Face, Nothing}(grid, boundary_conditions=boundary_conditions.v)
        velocities = (; u, v)
    end

    tracers = TracerFields(tracers, grid, boundary_conditions)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_salinity in a field 
    ice_salinity = field((Center, Center, Nothing), ice_salinity, grid)
    ice_density  = field((Center, Center, Nothing), ice_density, grid)

    # Construct prognostic fields if not provided
    ice_thickness     = Field{Center, Center, Nothing}(grid, boundary_conditions=boundary_conditions.h)
    ice_concentration = Field{Center, Center, Nothing}(grid, boundary_conditions=boundary_conditions.ℵ) 

    # Adding thickness and concentration if not there
    prognostic_fields = merge(tracers, (; h = ice_thickness, ℵ = ice_concentration))
    prognostic_fields = if ice_salinity isa ConstantField
        prognostic_fields 
    else
        merge(prognostic_fields, (; S = ice_salinity))
    end
    
    prognostic_fields = isnothing(dynamics) ? prognostic_fields : merge(prognostic_fields, velocities)

    # TODO: should we have ice thickness and concentration as part of the tracers or
    # just additional fields of the sea ice model?
    timestepper = ForwardEulerTimeStepper(grid, prognostic_fields)
    tracers     = merge(tracers, (; S = ice_salinity))

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

    arch = architecture(grid)

    return SeaIceModel(arch,
                       grid,
                       clock,
                       forcing, 
                       velocities,
                       tracers,
                       ice_thickness,
                       ice_concentration,
                       ice_density,
                       ice_consolidation_thickness,
                       ice_thermodynamics,
                       dynamics,
                       external_heat_fluxes,
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
                           fields(model.dynamics))

# TODO: make this correct
prognostic_fields(model::SIM) = merge((; h  = model.ice_thickness,
                                         ℵ  = model.ice_concentration),
                                      model.tracers,
                                      model.velocities)
