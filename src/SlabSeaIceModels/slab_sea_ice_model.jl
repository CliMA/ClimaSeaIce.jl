using ClimaSeaIce:
    PhaseTransitions,
    ForwardEulerTimestepper

using ClimaSeaIce.HeatBoundaryConditions:
    MeltingConstrainedFluxBalance,
    IceWaterThermalEquilibrium,
    PrescribedTemperature,
    FluxFunction,
    flux_summary

using Oceananigans.Utils: prettysummary
using Oceananigans.TimeSteppers: Clock, TimeStepper
using Oceananigans.Fields: field, Field, Center, ZeroField, ConstantField

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.Models: AbstractModel
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Simulations: reset!, initialize!, iteration
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Utils: prettytime

# TODO: move to Oceananigans
# import Oceananigans.Fields: field
# field(loc, a::Number, grid) = ConstantField(a)

struct SlabSeaIceModel{GR, CL, TS, IT, IC, ST, IS, U, D, UO, R, A, C, STF, SVF, TBC, CF, P, MIT} <: AbstractModel{TS}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # State
    thickness :: IT
    concentration :: IC
    top_surface_temperature :: ST
    salinity :: IS
    velocities :: U
    # Advection
    ocean_velocities :: UO
    rheology :: R
    advection :: A
    coriolis :: C
    # Boundary conditions
    external_heat_fluxes :: STF
    external_momentum_stress :: SVF
    heat_boundary_conditions :: TBC
    # Internal flux
    internal_heat_flux :: CF
    # Melting and freezing stuff
    phase_transitions :: P
    consolidation_thickness :: MIT
end

const SSIM = SlabSeaIceModel

Base.summary(model::SSIM) = "SlabSeaIceModel"
prettytime(model::SSIM) = prettytime(model.clock.time)
iteration(model::SSIM) = model.clock.iteration

function Base.show(io::IO, model::SSIM)
    grid = model.grid
    arch = architecture(grid)
    gridname = typeof(grid).name.wrapper
    timestr = string("(time = ", prettytime(model), ", iteration = ", iteration(model), ")")

    print(io, "SlabSeaIceModel{", typeof(arch), ", ", gridname, "}", timestr, '\n')
    print(io, "├── grid: ", summary(model.grid), '\n')
    print(io, "├── top_surface_temperature: ", summary(model.top_surface_temperature), '\n')
    print(io, "├── minimium_thickness: ", prettysummary(model.consolidation_thickness), '\n')
    print(io, "└── external_heat_fluxes: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_heat_fluxes.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_heat_fluxes.bottom, "     "))
    print(io, "└── external x-momentum stress: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_momentum_stress.u.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_momentum_stress.u.bottom, "     "))
    print(io, "└── external y-momentum stress: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_momentum_stress.v.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_momentum_stress.v.bottom, "     "))
end
         
reset!(::SSIM) = nothing
initialize!(::SSIM) = nothing
default_included_properties(::SSIM) = tuple(:grid)

fields(model::SSIM) = (h  = model.thickness,
                       ℵ  = model.concentration,
                       u  = model.velocities.u,
                       v  = model.velocities.v,
                       Tᵤ = model.top_surface_temperature,
                       Sᵢ = model.salinity)

# TODO: make this correct
prognostic_fields(model::SSIM) = fields(model)

"""
    SlabSeaIceModel(grid; kw...)

Pretty simple model for sea ice.
"""
function SlabSeaIceModel(grid;
                         clock                          = Clock{eltype(grid)}(0, 0, 1),
                         thickness                      = Field{Center, Center, Nothing}(grid),
                         consolidation_thickness        = 0.0, # m
                         concentration                  = Field{Center, Center, Nothing}(grid),
                         salinity                       = 0, # psu
                         tracers                        = (:h, :ℵ),
                         top_surface_temperature        = nothing,
                         top_heat_flux                  = nothing,
                         bottom_heat_flux               = 0,
                         top_u_stress                   = Field{Face, Center, Nothing}(grid),
                         top_v_stress                   = Field{Center, Face, Nothing}(grid),
                         bottom_u_stress                = Field{Face, Center, Nothing}(grid),
                         bottom_v_stress                = Field{Center, Face, Nothing}(grid),
                         velocities                     = nothing,
                         ocean_velocities               = (u = ZeroField(grid), v = ZeroField(grid)),
                         rheology                       = nothing,
                         advection                      = nothing,
                         coriolis                       = nothing,
                         top_heat_boundary_condition    = MeltingConstrainedFluxBalance(),
                         bottom_heat_boundary_condition = IceWaterThermalEquilibrium(),
                         internal_heat_flux             = ConductiveFlux(eltype(grid), conductivity=2),
                         phase_transitions              = PhaseTransitions(eltype(grid)))

    if isnothing(velocities) 
        u = XFaceField(grid)
        v = YFaceField(grid)
        velocities = (; u, v)
    end

    tracers = unique((tracers..., :h, :ℵ))

    if isnothing(salinity) # Prognostic salinity
        salinity = Field{Center, Center, Nothing}(grid)
        tracers = unique((tracers..., :S))
    else # Prescribed salinity (wrap salinity in a field)
        salinity = field((Center, Center, Nothing), salinity, grid)
    end

    # Only one time-stepper is supported currently
    timestepper = TimeStepper(:QuasiAdamsBashforth2, grid, tracernames(tracers);
                              Gⁿ = SlabSeaIceModelTendencyFields(grid, tracernames(:h, :ℵ)),
                              G⁻ = SlabSeaIceModelTendencyFields(grid, tracernames(:h, :ℵ)))

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    consolidation_thickness = field((Center, Center, Nothing), consolidation_thickness, grid)

    # Construct an internal heat flux function that captures the liquidus and
    # bottom boundary condition.
    parameters = (flux = internal_heat_flux,
                  liquidus = phase_transitions.liquidus,
                  bottom_heat_boundary_condition = bottom_heat_boundary_condition)

    internal_heat_flux_function = FluxFunction(slab_internal_heat_flux;
                                                  parameters,
                                                  top_temperature_dependent=true)

    # Construct default top heat flux if one is not provided
    if isnothing(top_heat_flux)
        if top_heat_boundary_condition isa PrescribedTemperature  
            # Default: external top flux is in equilibrium with internal fluxes
            top_heat_flux = internal_heat_flux_function
        else
            # Default: no external top surface flux
            top_heat_flux = 0
        end
    end

    # Construct default top temperature if one is not provided
    if isnothing(top_surface_temperature)
        # Check top boundary condition
        if top_heat_boundary_condition isa PrescribedTemperature  
            top_surface_temperature = top_heat_boundary_condition.temperature 
        else # build the default
            top_surface_temperature = Field{Center, Center, Nothing}(grid)
        end
    end

    # Convert to `field` (does nothing if it's already a Field)
    top_surface_temperature = field((Center, Center, Nothing), top_surface_temperature, grid)

    # Package the external fluxes and boundary conditions
    external_heat_fluxes = (top = top_heat_flux,    
                         bottom = bottom_heat_flux) 

    external_momentum_stress = (; u = (; top = top_u_stress,
                                         bottom = bottom_u_stress),
                                  v = (; top = top_v_stress, 
                                         bottom = bottom_v_stress))

    heat_boundary_conditions = (top = top_heat_boundary_condition,
                             bottom = bottom_heat_boundary_condition)

    return SlabSeaIceModel(grid,
                           clock,
                           timestepper,
                           thickness,
                           concentration,
                           top_surface_temperature,
                           salinity,
                           velocities,
                           ocean_velocities,
                           rheology,
                           advection,
                           coriolis,
                           external_heat_fluxes,
                           external_momentum_stress,
                           heat_boundary_conditions,
                           internal_heat_flux_function,
                           phase_transitions,
                           consolidation_thickness)
end

function set!(model::SSIM; h=nothing, ℵ=nothing)
    !isnothing(h) && set!(model.thickness, h)
    !isnothing(ℵ) && set!(model.concentration, ℵ)
    return nothing
end

function SlabSeaIceModelTendencyFields(grid, tracer_names)
    u = XFaceField(grid)
    v = YFaceField(grid)
    tracers = TracerFields(tracer_names, grid)
    return merge((u=u, v=v), tracers)
end