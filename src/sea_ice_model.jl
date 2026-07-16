using Oceananigans: tupleit, tracernames
using Oceananigans.Advection: materialize_advection
using Oceananigans.BoundaryConditions: FieldBoundaryConditions, BoundaryCondition, Zipper,
                                       regularize_field_boundary_conditions
using Oceananigans.Fields: TracerFields, ConstantField
using Oceananigans.Forcings: model_forcing
using Oceananigans.Models: update_model_field_time_series!
using Oceananigans.OutputReaders: FieldTimeSeries
using Oceananigans.Simulations: iteration
using Oceananigans.TimeSteppers: TimeStepper
using Oceananigans.Utils: prettysummary

using .SeaIceDynamics: materialize_solver, maybe_extended_grid
using .SeaIceThermodynamics: PrescribedTemperature, FluxFunction,
                             PhaseTransitions, internal_flux_function,
                             writable_top_surface_temperature
using .SeaIceThermodynamics.HeatBoundaryConditions: flux_summary

@inline instantiate(T::DataType) = T()
@inline instantiate(T) = T

struct SeaIceModel{GR, TD, SNT, D, TS, CL, U, T, IT, IC, SNH, ID, SND, PT, CT, SP, MFX, STF, A, F, Arch} <: AbstractModel{TS, Arch}
    architecture :: Arch
    grid :: GR
    clock :: CL
    forcing :: F
    # Prognostic State
    velocities :: U
    tracers :: T
    ice_thickness :: IT
    ice_concentration :: IC
    snow_thickness :: SNH
    sea_ice_density :: ID
    snow_density :: SND
    ice_consolidation_thickness :: CT
    # Shared thermodynamic parameters
    phase_transitions :: PT
    # Thermodynamics
    ice_thermodynamics :: TD
    snow_thermodynamics :: SNT
    # Dynamics
    dynamics :: D
    # External boundary conditions
    external_heat_fluxes :: STF
    snowfall :: SP
    # Diagnostics (kg m⁻² s⁻¹)
    mass_fluxes :: MFX
    # Numerics
    timestepper :: TS
    advection :: A
end

assumed_sea_ice_field_location(name) = name === :u  ? (Face,   Center, Nothing) :
                                       name === :v  ? (Center, Face,   Nothing) :
                                                      (Center, Center, Nothing)

function default_sea_ice_boundary_conditions(grid, name)
    bcs = FieldBoundaryConditions(grid, instantiate.(assumed_sea_ice_field_location(name)))
    if (name === :u || name === :v) && bcs.north isa BoundaryCondition && bcs.north.classification isa Zipper
        north = BoundaryCondition(bcs.north.classification, - bcs.north.condition)
        bcs = FieldBoundaryConditions(bcs.west, bcs.east, bcs.south, north, bcs.bottom, bcs.top, bcs.immersed)
    end
    return bcs
end

"""
    SeaIceModel(grid;
                clock                       = Clock{eltype(grid)}(time = 0),
                ice_consolidation_thickness = 0.05, # m
                ice_salinity                = 0,    # psu
                sea_ice_density             = 900,  # kg m⁻³, bulk sea-ice
                snow_density                = 330,  # kg m⁻³, bulk snow
                phase_transitions           = PhaseTransitions(eltype(grid)),
                top_heat_flux               = nothing,
                bottom_heat_flux            = 0,    # W m⁻²
                velocities                  = nothing,
                advection                   = nothing,
                tracers                     = (),
                timestepper                 = :SplitRungeKutta3,
                boundary_conditions         = NamedTuple(),
                ice_thermodynamics          = sea_ice_slab_thermodynamics(grid),
                snow_thermodynamics         = nothing,
                snowfall                    = 0,
                dynamics                    = nothing,
                forcing                     = NamedTuple())

Construct a sea-ice model on `grid` with slab thermodynamics, optional momentum
equations, and optional snow thermodynamics.

The model evolves sea-ice thickness `h`, concentration `ℵ`, optional salinity
`S`, optional snow thickness `hs`, and, when `dynamics` is provided, the
horizontal ice velocity components `u` and `v`.

The field `mass_fluxes` holds three diagnostics (kg m⁻² s⁻¹) written by the thermodynamic step:
`mass_fluxes.thermodynamics.ice` and `mass_fluxes.thermodynamics.snow`, the per-step mass tendencies
of the ice and snow slabs (their sum is the mass exchanged with the ocean), and
`mass_fluxes.intercepted_snowfall`, the snowfall absorbed by the ice. These are excluded from
`fields(model)` and the prognostic state.

Arguments
=========

- `grid`: The computational grid.

Keyword Arguments
=================

- `clock`: Model clock. Defaults to `Clock{eltype(grid)}(time = 0)`.
- `ice_consolidation_thickness`: Threshold thickness above which the slab is
                                 treated as consolidated ice. Default: 0.05 meters.
- `ice_salinity`: Sea-ice salinity field or constant. When non-constant it is
                  included in the prognostic state as `S`. Default: 0 psu.
- `sea_ice_density`: Bulk sea-ice density. Default: 900 kg m⁻³.
- `snow_density`: Bulk snow density. Default: 330 kg m⁻³.
- `phase_transitions`: Thermodynamic phase-transition parameters shared by the
                       ice and snow slabs. Default:
                       `PhaseTransitions(eltype(grid))`.
- `top_heat_flux`: External top heat flux, or tuple of fluxes, passed to the
                   ice upper boundary condition. Default: `nothing`.
- `bottom_heat_flux`: External bottom heat flux passed to the ice lower
                      boundary condition. Default: `0` W m⁻².
- `velocities`: Pre-existing velocity fields. If omitted, they are allocated by
                the constructor, possibly on an extended grid required by the
                selected dynamics solver. Default: `nothing`.
- `advection`: Advection scheme for prognostic tracers. Default: `nothing`.
- `tracers`: Additional prognostic tracers.
- `timestepper`: Time stepper specification passed to `TimeStepper`. Default:
                 `:SplitRungeKutta3`.
- `boundary_conditions`: Boundary conditions for allocated prognostic fields.
- `ice_thermodynamics`: Ice thermodynamics model. Default:
                        `sea_ice_slab_thermodynamics(grid)`.
- `snow_thermodynamics`: Optional snow thermodynamics model. When provided, the
                         model allocates prognostic snow thickness `hs`.
- `snowfall`: Snowfall forcing applied only when `snow_thermodynamics` is
              present. May be a constant, `Field`, or `FieldTimeSeries`.
- `dynamics`: Optional sea-ice dynamics model, for example,
              `SeaIceMomentumEquation(grid)`.
- `forcing`: Additional model forcing passed through `model_forcing`.
"""
function SeaIceModel(grid;
                     clock                       = Clock{eltype(grid)}(time = 0),
                     ice_consolidation_thickness = 0.05, # m
                     ice_salinity                = 0,    # psu
                     sea_ice_density             = 900,  # kg m⁻³, bulk sea-ice
                     snow_density                = 330,  # kg m⁻³, bulk snow
                     phase_transitions           = PhaseTransitions(eltype(grid)),
                     top_heat_flux               = nothing,
                     bottom_heat_flux            = 0,    # W m⁻²
                     velocities                  = nothing,
                     advection                   = nothing,
                     tracers                     = (),
                     timestepper                 = :SplitRungeKutta3,
                     boundary_conditions         = NamedTuple(),
                     ice_thermodynamics          = sea_ice_slab_thermodynamics(grid),
                     snow_thermodynamics         = nothing,
                     snowfall                    = 0,
                     dynamics                    = nothing,
                     forcing                     = NamedTuple())

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_consolidation_thickness in a field
    ice_consolidation_thickness = field((Center, Center, Nothing), ice_consolidation_thickness, grid)

    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)

    # Next, we form a list of default boundary conditions:
    field_names = (:u, :v, :h, :ℵ, :S, tracernames(tracers)...)

    bc_tuple = Tuple(default_sea_ice_boundary_conditions(grid, name) for name in field_names)
    default_boundary_conditions = NamedTuple{field_names}(bc_tuple)

    # Then we merge specified, embedded, and default boundary conditions. Specified boundary conditions
    # have precedence, followed by embedded, followed by default.
    boundary_conditions = merge(default_boundary_conditions, boundary_conditions)
    boundary_conditions = regularize_field_boundary_conditions(boundary_conditions, grid, field_names)

    if isnothing(velocities)
        # Extend the halos for the velocity fields if the dynamics is
        # an split explicit momentum equation on a Distributed grid
        velocity_grid = maybe_extended_grid(dynamics, grid)
        dynamics = materialize_solver(dynamics, velocity_grid)
        u = Field{Face, Center, Nothing}(velocity_grid, boundary_conditions=boundary_conditions.u)
        v = Field{Center, Face, Nothing}(velocity_grid, boundary_conditions=boundary_conditions.v)
        velocities = (; u, v)
    else
        velocity_grid = velocities.u.grid
    end

    tracers = TracerFields(tracers, grid, boundary_conditions)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    ice_salinity    = field((Center, Center, Nothing), ice_salinity,    grid)
    sea_ice_density = field((Center, Center, Nothing), sea_ice_density, grid)
    snow_density    = field((Center, Center, Nothing), snow_density,    grid)

    # Thickness and concentration need to be on the velocity_grid because rheology needs both `h` and `ℵ`
    # _inside_ the halos when running in an extended distributed grid. In serial cases and for solvers
    # other than split explicit velocity_grid == grid.
    ice_thickness     = Field{Center, Center, Nothing}(velocity_grid, boundary_conditions=boundary_conditions.h)
    ice_concentration = Field{Center, Center, Nothing}(velocity_grid, boundary_conditions=boundary_conditions.ℵ)

    # Snow thickness
    snow_thickness = if isnothing(snow_thermodynamics)
        nothing
    else
        Field{Center, Center, Nothing}(grid)
    end

    # Wrap snowfall in a field (but preserve FieldTimeSeries)
    if !(snowfall isa FieldTimeSeries)
        snowfall = field((Center, Center, Nothing), snowfall, grid)
    end

    # Adding sea-ice thickness and concentration to prognostic fields if not there already
    prognostic_fields = merge(tracers, (; h = ice_thickness, ℵ = ice_concentration))

    # Add snow thickness to prognostic fields when present
    prognostic_fields = if !isnothing(snow_thickness)
        merge(prognostic_fields, (; hs = snow_thickness))
    else
        prognostic_fields
    end

    prognostic_fields = if ice_salinity isa ConstantField
        prognostic_fields
    else
        merge(prognostic_fields, (; S = ice_salinity))
    end

    prognostic_fields = isnothing(dynamics) ? prognostic_fields : merge(prognostic_fields, velocities)

    # TODO: should we have ice thickness and concentration as part of the tracers or
    # just additional fields of the sea ice model?
    tracers = merge(tracers, (; S = ice_salinity))
    timestepper = TimeStepper(timestepper, grid, prognostic_fields)

    # The layered (snow + ice) step writes the ice top surface temperature, so it
    # must be writable when snow is present; bare-ice models keep their field as-is.
    if !isnothing(ice_thermodynamics) && !isnothing(snow_thermodynamics)
        ice_thermodynamics = writable_top_surface_temperature(ice_thermodynamics, grid)
    end

    if !isnothing(ice_thermodynamics)
        if isnothing(top_heat_flux)
            if isnothing(snow_thermodynamics) &&
               ice_thermodynamics.heat_boundary_conditions.top isa PrescribedTemperature
                # Default: external top flux is in equilibrium with internal fluxes.
                # Build a FluxFunction wrapper using the model's shared liquidus.
                top_heat_flux = internal_flux_function(ice_thermodynamics.internal_heat_flux,
                                                       phase_transitions.liquidus,
                                                       ice_thermodynamics.heat_boundary_conditions.bottom)
            else
                # Default: no external top surface flux
                top_heat_flux = 0
            end
        end
    end

    model_fields = isnothing(dynamics) ? prognostic_fields : merge(prognostic_fields, fields(dynamics))
    forcing = model_forcing(forcing, model_fields, prognostic_fields)

    # Fill any settings in advection scheme that might have been deferred until
    # the grid and backend is known
    advection = materialize_advection(advection, grid)

    # Package the external fluxes and boundary conditions
    external_heat_fluxes = (top = top_heat_flux,
                            bottom = bottom_heat_flux)

    mass_fluxes = (thermodynamics = (ice  = Field{Center, Center, Nothing}(grid),
                                     snow = Field{Center, Center, Nothing}(grid)),
                   intercepted_snowfall = Field{Center, Center, Nothing}(grid))

    arch = architecture(grid)

    return SeaIceModel(arch,
                       grid,
                       clock,
                       forcing,
                       velocities,
                       tracers,
                       ice_thickness,
                       ice_concentration,
                       snow_thickness,
                       sea_ice_density,
                       snow_density,
                       ice_consolidation_thickness,
                       phase_transitions,
                       ice_thermodynamics,
                       snow_thermodynamics,
                       dynamics,
                       external_heat_fluxes,
                       snowfall,
                       mass_fluxes,
                       timestepper,
                       advection)
end

const SIM = SeaIceModel

function Oceananigans.Fields.set!(model::SIM; h=nothing, ℵ=nothing, hs=nothing, u=nothing, v=nothing)
    !isnothing(h)  && set!(model.ice_thickness, h)
    !isnothing(ℵ)  && set!(model.ice_concentration, ℵ)
    !isnothing(u)  && set!(model.velocities.u, u)
    !isnothing(v)  && set!(model.velocities.v, v)

    if !isnothing(hs)
        if isnothing(model.snow_thickness)
            throw(ArgumentError("Cannot set snow thickness `hs` on a SeaIceModel without snow (model.snow_thickness is nothing)."))
        end
        set!(model.snow_thickness, hs)
    end

    return nothing
end

Oceananigans.Fields.set!(model::SIM, new_clock::Clock) = set!(model.clock, new_clock)

Oceananigans.Architectures.architecture(model::SIM) = model.architecture
Oceananigans.Simulations.iteration(model::SIM) = model.clock.iteration
Oceananigans.Utils.prettytime(model::SIM) = prettytime(model.clock.time)

function Base.summary(model::SIM)
    A = Base.summary(architecture(model))
    G = nameof(typeof(model.grid))
    return string("SeaIceModel{$A, $G}",
                  "(time = ", prettytime(model.clock.time),
                  ", iteration = ", prettysummary(model.clock.iteration), ")")
end

function Base.show(io::IO, model::SIM)
    grid = model.grid
    arch = architecture(grid)
    gridname = typeof(grid).name.wrapper
    timestr = string("(time = ", prettytime(model), ", iteration = ", iteration(model), ")")

    print(io, "SeaIceModel{", typeof(arch), ", ", gridname, "}", timestr, '\n')
    print(io, "├── grid: ", summary(model.grid), '\n')
    print(io, "├── timestepper: ", summary(model.timestepper), '\n')
    print(io, "├── ice_thermodynamics: ", summary(model.ice_thermodynamics), '\n')
    print(io, "├── snow_thermodynamics: ", summary(model.snow_thermodynamics), '\n')
    print(io, "├── advection: ", summary(model.advection), '\n')
    print(io, "└── external_heat_fluxes: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_heat_fluxes.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_heat_fluxes.bottom, "     "))
end

Oceananigans.initialize!(::SIM) = nothing
Oceananigans.OutputWriters.default_included_properties(::SIM) = [:grid]
Oceananigans.OutputWriters.checkpointer_address(::SIM) = "SeaIceModel"
Oceananigans.TimeSteppers.reset!(::SIM) = nothing

snow_fields(::Nothing) = NamedTuple()
snow_fields(hs) = (; hs)

component_fields(component) = fields(component)
component_fields(::Nothing) = NamedTuple()

component_prognostic_fields(component) = prognostic_fields(component)
component_prognostic_fields(::Nothing) = NamedTuple()
component_prognostic_fields(model, component) = prognostic_fields(model, component)
component_prognostic_fields(model, ::Nothing) = NamedTuple()

Oceananigans.fields(model::SIM) = merge((; h  = model.ice_thickness,
                                           ℵ  = model.ice_concentration,
                                           ρi = model.sea_ice_density),
                                        snow_fields(model.snow_thickness),
                                        model.tracers,
                                        model.velocities,
                                        component_fields(model.ice_thermodynamics),
                                        component_fields(model.dynamics))

Oceananigans.prognostic_fields(model::SIM) = merge((; h  = model.ice_thickness,
                                                      ℵ  = model.ice_concentration),
                                                   snow_fields(model.snow_thickness),
                                                   component_prognostic_fields(model, model.dynamics),
                                                   component_prognostic_fields(model.ice_thermodynamics))

function Oceananigans.TimeSteppers.update_state!(model::SIM, callbacks=[])

    foreach(prognostic_fields(model)) do field
        mask_immersed_field_xy!(field, k=size(model.grid, 3))
        fill_halo_regions!(field, model.clock, fields(model))
    end

    Nz = size(model.grid, 3)
    mask_immersed_field_xy!(model.mass_fluxes.thermodynamics.ice,   k=Nz)
    mask_immersed_field_xy!(model.mass_fluxes.thermodynamics.snow,  k=Nz)
    mask_immersed_field_xy!(model.mass_fluxes.intercepted_snowfall, k=Nz)

    update_model_field_time_series!(model, model.clock)

    return nothing
end

function Oceananigans.Models.update_model_field_time_series!(model::SIM, clock::Clock)
    time = Time(clock.time)

    possible_fts = (model.tracers, model.external_heat_fluxes, model.snowfall, model.dynamics)
    time_series_tuple = extract_field_time_series(possible_fts)
    time_series_tuple = flattened_unique_values(time_series_tuple)

    for fts in time_series_tuple
        update_field_time_series!(fts, time)
    end

    return nothing
end

#####
##### Checkpointing
#####

function Oceananigans.prognostic_state(model::SIM)
    return (clock = prognostic_state(model.clock),
            velocities = prognostic_state(model.velocities),
            ice_thickness = prognostic_state(model.ice_thickness),
            ice_concentration = prognostic_state(model.ice_concentration),
            snow_thickness = prognostic_state(model.snow_thickness),
            tracers = prognostic_state(model.tracers),
            timestepper = prognostic_state(model.timestepper),
            ice_thermodynamics = prognostic_state(model.ice_thermodynamics),
            snow_thermodynamics = prognostic_state(model.snow_thermodynamics),
            dynamics = prognostic_state(model.dynamics),
            mass_fluxes = prognostic_state(model.mass_fluxes))
end

function Oceananigans.restore_prognostic_state!(model::SIM, state)
    restore_prognostic_state!(model.clock, state.clock)
    restore_prognostic_state!(model.velocities, state.velocities)
    restore_prognostic_state!(model.ice_thickness, state.ice_thickness)
    restore_prognostic_state!(model.ice_concentration, state.ice_concentration)
    restore_prognostic_state!(model.snow_thickness, state.snow_thickness)
    restore_prognostic_state!(model.tracers, state.tracers)
    restore_prognostic_state!(model.timestepper, state.timestepper)
    restore_prognostic_state!(model.ice_thermodynamics, state.ice_thermodynamics)
    restore_prognostic_state!(model.snow_thermodynamics, state.snow_thermodynamics)
    restore_prognostic_state!(model.dynamics, state.dynamics)
    if haskey(state, :mass_fluxes) # Backwards compatible (TODO: remove)
        restore_prognostic_state!(model.mass_fluxes, state.mass_fluxes)
    end
    return model
end

Oceananigans.restore_prognostic_state!(::SIM, ::Nothing) = nothing
