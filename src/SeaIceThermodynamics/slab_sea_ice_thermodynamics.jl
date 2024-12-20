import Oceananigans: fields

struct SlabSeaIceThermodynamics{ST, HBC, CF, P, MIT}
    top_surface_temperature :: ST
    heat_boundary_conditions :: HBC
    # Internal flux
    internal_heat_flux :: CF
    # Melting and freezing stuff
    phase_transitions :: P
    ice_consolidation_thickness :: MIT
end

Adapt.adapt_structure(to, t::SlabSeaIceThermodynamics) = 
    SlabSeaIceThermodynamics(Adapt.adapt(to, t.top_surface_temperature),
                             Adapt.adapt(to, t.heat_boundary_conditions),
                             Adapt.adapt(to, t.internal_heat_flux),
                             Adapt.adapt(to, t.phase_transitions),
                             Adapt.adapt(to, t.ice_consolidation_thickness))

const SSIT = SlabSeaIceThermodynamics

Base.summary(therm::SSIT) = "SlabThermodynamics"

function Base.show(io::IO, therm::SSIT)
    print(io, "SlabSeaIceThermodynamics", '\n')
    print(io, "├── top_surface_temperature: ", summary(therm.top_surface_temperature), '\n')
    print(io, "└── minimium_ice_thickness: ", prettysummary(therm.ice_consolidation_thickness), '\n')
end
       
fields(therm::SSIT) = (; Tu = therm.top_surface_temperature)

"""
    SlabSeaIceThermodynamics(grid; kw...)

Pretty simple model for sea ice.
"""
function SlabSeaIceThermodynamics(grid;
                                  ice_consolidation_thickness    = 0.0, # m
                                  top_surface_temperature        = nothing,
                                  top_heat_boundary_condition    = MeltingConstrainedFluxBalance(),
                                  bottom_heat_boundary_condition = IceWaterThermalEquilibrium(),
                                  # Default internal flux: thermal conductivity of 2 kg m s⁻³ K⁻¹, appropriate for freshwater ice
                                  internal_heat_flux             = ConductiveFlux(eltype(grid), conductivity=2),
                                  phase_transitions              = PhaseTransitions(eltype(grid)))

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_consolidation_thickness in a field
    ice_consolidation_thickness = field((Center, Center, Nothing), ice_consolidation_thickness, grid)

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

    return SlabSeaIceThermodynamics(top_surface_temperature,
                                    heat_boundary_conditions,
                                    internal_heat_flux_function,
                                    phase_transitions,
                                    ice_consolidation_thickness)
end

