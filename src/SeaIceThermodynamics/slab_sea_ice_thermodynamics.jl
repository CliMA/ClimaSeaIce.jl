import Oceananigans: fields

struct ProportionalEvolution end

struct SlabSeaIceThermodynamics{ST, HBC, CF, GT, P, CE}
    top_surface_temperature :: ST
    heat_boundary_conditions :: HBC
    # Internal flux
    internal_heat_flux :: CF
    thermodynamic_tendency :: GT
    # Melting and freezing stuff
    phase_transitions :: P
    # Rules to evolve concentration
    concentration_evolution :: CE
end

Adapt.adapt_structure(to, t::SlabSeaIceThermodynamics) = 
    SlabSeaIceThermodynamics(Adapt.adapt(to, t.top_surface_temperature),
                             Adapt.adapt(to, t.heat_boundary_conditions),
                             Adapt.adapt(to, t.internal_heat_flux),
                             Adapt.adapt(to, t.thermodynamic_tendency),
                             Adapt.adapt(to, t.phase_transitions),
                             Adapt.adapt(to, t.concentration_evolution))

const SSIT = SlabSeaIceThermodynamics

Base.summary(therm::SSIT) = "SlabThermodynamics"

function Base.show(io::IO, therm::SSIT)
    print(io, "SlabSeaIceThermodynamics", '\n')
    print(io, "└── top_surface_temperature: ", summary(therm.top_surface_temperature))
end
       
fields(therm::SSIT) = (; Tu = therm.top_surface_temperature)

"""
    SlabSeaIceThermodynamics(grid; kw...)

Pretty simple model for sea ice.
"""
function SlabSeaIceThermodynamics(grid;
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

    thermodynamic_tendency = Field{Center, Center, Nothing}(grid)
    heat_boundary_conditions = (top = top_heat_boundary_condition,
                                bottom = bottom_heat_boundary_condition)

    return SlabSeaIceThermodynamics(top_surface_temperature,
                                    heat_boundary_conditions,
                                    internal_heat_flux_function,
                                    thermodynamic_tendency,
                                    phase_transitions,
                                    concentration_evolution)
end

