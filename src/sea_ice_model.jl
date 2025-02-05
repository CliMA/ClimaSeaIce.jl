using Oceananigans.Fields: TracerFields
using Oceananigans.TimeSteppers: TimeStepper
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions
using ClimaSeaIce.SeaIceThermodynamics: PrescribedTemperature
using Oceananigans: tupleit, tracernames
using Oceananigans.BoundaryConditions: regularize_field_boundary_conditions
using Oceananigans.Forcings: model_forcing
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: flux_summary
using ClimaSeaIce.Rheologies: required_prognostic_tracers

struct SeaIceModel{GR, TD, D, TS, CL, U, T, IT, IC, ID, CT, STF, A, F} <: AbstractModel{TS}
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

function SeaIceModel(grid;
                     clock                       = Clock{eltype(grid)}(time = 0),
                     ice_thickness               = nothing,
                     ice_consolidation_thickness = 0.05, # m
                     ice_concentration           = nothing,
                     ice_salinity                = 0, # psu
                     ice_density                 = 900, # kg m⁻³
                     top_heat_flux               = nothing,
                     bottom_heat_flux            = 0,
                     velocities                  = nothing,
                     timestepper                 = :QuasiAdamsBashforth2,
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

    # Next, we form a list of default boundary conditions:
    field_names = (:u, :v, :h, :ℵ, :S, tracernames(tracers)...)
    default_boundary_conditions = NamedTuple{field_names}(Tuple(FieldBoundaryConditions()
                                                          for name in field_names))

    # Then we merge specified, embedded, and default boundary conditions. Specified boundary conditions
    # have precedence, followed by embedded, followed by default.
    boundary_conditions = merge(default_boundary_conditions, boundary_conditions)
    boundary_conditions = regularize_field_boundary_conditions(boundary_conditions, grid, field_names)
    
    if isnothing(velocities) 
        u = Field{Face, Face, Center}(grid) #, boundary_conditions=boundary_conditions.u)
        v = Field{Face, Face, Center}(grid) #, boundary_conditions=boundary_conditions.v)
        velocities = (; u, v)
    end

    tracers = TracerFields(tracers, grid, boundary_conditions)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_salinity in a field 
    ice_salinity = field((Center, Center, Center), ice_salinity, grid)
    ice_density  = field((Center, Center, Center), ice_density,  grid)

    # Construct prognostic fields if not provided
    ice_thickness = isnothing(ice_thickness) ?  Field{Center, Center, Center}(grid, boundary_conditions=boundary_conditions.h) : ice_thickness
    ice_concentration = isnothing(ice_concentration) ? Field{Center, Center, Center}(grid, boundary_conditions=boundary_conditions.ℵ) : ice_concentration

    rheology_tracers = required_prognostic_tracers(dynamics.rheology, grid)

    # Adding thickness and concentration if not there
    prognostic_fields = merge(tracers, rheology_tracers, (; h = ice_thickness, ℵ = ice_concentration))
    prognostic_fields = if ice_salinity isa ConstantField 
        prognostic_fields 
    else
        merge(prognostic_fields, (; S = ice_salinity))
    end
    
    prognostic_fields = isnothing(dynamics) ? prognostic_fields : merge(prognostic_fields, velocities)

    # TODO: should we have ice thickness and concentration as part of the tracers or
    # just additional fields of the sea ice model?
    tracers = merge(tracers, rheology_tracers, (; S = ice_salinity))
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

    return SeaIceModel(grid,
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

@kernel function _set_minium_ice_thickness!(h, ℵ, hmin)
    i, j = @index(Global, NTuple)

    @inbounds begin
        h⁺ = h[i, j, 1]
        ℵ⁺ = ℵ[i, j, 1]
        h⁻ = hmin[i, j, 1]

        ht, ℵt = cap_ice_thickness(h⁺, h⁻, ℵ⁺)

        ℵ[i, j, 1] = ℵt
        h[i, j, 1] = ht
    end
end

function set!(model::SIM; h=nothing, ℵ=nothing, kwargs...)
    grid = model.grid
    arch = architecture(model)

    !isnothing(h) && set!(model.ice_thickness, h)
    !isnothing(ℵ) && set!(model.ice_concentration, ℵ)

        #We cap the ice to the consolidation thickness
    launch!(arch, grid, :xy, _set_minium_ice_thickness!, 
            model.ice_thickness,
            model.ice_concentration,
            model.ice_consolidation_thickness)

    # All the other tracers
    for (fldname, value) in kwargs
        if fldname ∈ propertynames(model.velocities)
            ϕ = getproperty(model.velocities, fldname)
        elseif fldname ∈ propertynames(model.tracers)
            ϕ = getproperty(model.tracers, fldname)
        else
            throw(ArgumentError("name $fldname not found in model.velocities, model.tracers, or model.free_surface"))
        end

        set!(ϕ, value)
    end

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
