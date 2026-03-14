import Oceananigans: fields, prognostic_fields, prognostic_state, restore_prognostic_state!

struct ProportionalEvolution end

struct SlabThermodynamics{ST, HBC, CF, P, CE}
    top_surface_temperature :: ST
    heat_boundary_conditions :: HBC
    internal_heat_flux :: CF
    phase_transitions :: P
    concentration_evolution :: CE
end

Adapt.adapt_structure(to, t::SlabThermodynamics) =
    SlabThermodynamics(Adapt.adapt(to, t.top_surface_temperature),
                       Adapt.adapt(to, t.heat_boundary_conditions),
                       Adapt.adapt(to, t.internal_heat_flux),
                       Adapt.adapt(to, t.phase_transitions),
                       Adapt.adapt(to, t.concentration_evolution))

const SSIT = SlabThermodynamics

# Backward compatibility alias
const SlabSeaIceThermodynamics = SlabThermodynamics

"""
    SlabSnowThermodynamics(grid; kw...)

Construct a `SlabThermodynamics` with default parameters appropriate for snow:
conductivity = 0.31 W/(m K), density = 330 kg/m³, heat capacity = 2090 J/(kg K),
and latent heat = 334000 J/kg.
"""
function SlabSnowThermodynamics(grid;
                                conductivity          = 0.31,
                                density               = 330,
                                heat_capacity         = 2090,
                                reference_latent_heat = 334e3,
                                kw...)

    FT = eltype(grid)
    internal_heat_flux = ConductiveFlux(FT, conductivity = conductivity)
    phase_transitions  = PhaseTransitions(FT;
                                          density               = density,
                                          heat_capacity         = heat_capacity,
                                          reference_latent_heat = reference_latent_heat)

    return SlabThermodynamics(grid; internal_heat_flux, phase_transitions, kw...)
end

Base.summary(therm::SSIT) = "SlabThermodynamics"

function Base.show(io::IO, therm::SSIT)
    print(io, "SlabThermodynamics", '\n')
    print(io, "└── top_surface_temperature: ", summary(therm.top_surface_temperature))
end

fields(therm::SSIT) = (; Tu = therm.top_surface_temperature)
prognostic_fields(therm::SSIT) = NamedTuple()

"""
    SlabThermodynamics(grid; kw...)

Pretty simple model for sea ice.
"""
function SlabThermodynamics(grid;
                            top_surface_temperature        = nothing,
                            top_heat_boundary_condition    = MeltingConstrainedFluxBalance(),
                            bottom_heat_boundary_condition = IceWaterThermalEquilibrium(),
                            # Default internal flux: thermal conductivity of 2 kg m s⁻³ K⁻¹, appropriate for freshwater ice
                            internal_heat_flux             = ConductiveFlux(eltype(grid), conductivity=2),
                            phase_transitions              = PhaseTransitions(eltype(grid)),
                            concentration_evolution        = ProportionalEvolution())

    # Construct an internal heat flux function that captures the liquidus and
    # bottom boundary condition.
    parameters = (flux = internal_heat_flux,
                  liquidus = phase_transitions.liquidus,
                  bottom_heat_boundary_condition = bottom_heat_boundary_condition)

    internal_heat_flux_function = FluxFunction(slab_internal_heat_flux;
                                               parameters,
                                               top_temperature_dependent=true)

    if top_heat_boundary_condition isa PrescribedTemperature
        if !isnothing(top_surface_temperature)
            msg = "You cannot provide a redundant top_surface_temperature when using \
                   PrescribedTemperature top_heat_boundary_condition."
            throw(ArgumentError(msg))
        else
            # Convert to `field` (does nothing if it's already a Field)
            top_surface_temperature = top_heat_boundary_condition.temperature
            top_surface_temperature = field((Center, Center, Nothing), top_surface_temperature, grid)
        end
    else
        top_surface_temperature = Field{Center, Center, Nothing}(grid)
    end

    heat_boundary_conditions = (top = top_heat_boundary_condition,
                                bottom = bottom_heat_boundary_condition)

    return SlabThermodynamics(top_surface_temperature,
                              heat_boundary_conditions,
                              internal_heat_flux_function,
                              phase_transitions,
                              concentration_evolution)
end

#####
##### Checkpointing
#####

function prognostic_state(therm::SlabThermodynamics)
    return (top_surface_temperature = prognostic_state(therm.top_surface_temperature),)
end

function restore_prognostic_state!(therm::SlabThermodynamics, state)
    restore_prognostic_state!(therm.top_surface_temperature, state.top_surface_temperature)
    return therm
end
