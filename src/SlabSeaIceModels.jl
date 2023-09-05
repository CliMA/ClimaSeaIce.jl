module SlabSeaIceModels

using ClimaSeaIce:
    PhaseTransitions,
    latent_heat,
    ForwardEulerTimestepper,
    melting_temperature

using ClimaSeaIce.ThermalBoundaryConditions:
    MeltingConstrainedFluxBalance,
    IceWaterThermalEquilibrium,
    PrescribedTemperature,
    FluxFunction,
    SurfaceTemperatureDependent,
    bottom_temperature,
    top_temperature,
    top_flux_imbalance,
    bottom_flux_imbalance,
    flux_summary

# using RootSolvers: find_zero

using Oceananigans.Architectures: architecture
using Oceananigans.Fields: Field, Center, ZeroField, ConstantField
using Oceananigans.TimeSteppers: Clock, tick!

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: field, set!
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!, initialize!
import Oceananigans.OutputWriters: default_included_properties
import Oceananigans.Utils: prettytime
import Oceananigans.Simulations: iteration

# TODO: move to Oceananigans
field(loc, a::Number, grid) = ConstantField(a)

struct SlabSeaIceModel{GR, CL, TS, IT, ST, IS, STF, TBC, CF, P}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # State
    ice_thickness :: IT
    top_temperature :: ST
    ice_salinity :: IS
    # Boundary conditions
    external_thermal_fluxes :: STF
    thermal_boundary_conditions :: TBC
    # Internal flux
    internal_thermal_flux :: CF
    # Melting and freezing stuff
    phase_transitions :: P
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
    print(io, "├── top_temperature: ", summary(model.top_temperature), '\n')
    print(io, "└── external_thermal_fluxes: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_thermal_fluxes.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_thermal_fluxes.bottom, "     "))
end
         
initialize!(::SSIM) = nothing
default_included_properties(::SSIM) = tuple(:grid)

fields(model::SSIM) = (h = model.ice_thickness,
                       Tᵤ = model.top_temperature,
                       Sᵢ = model.ice_salinity)

# TODO: make this correct
prognostic_fields(model::SSIM) = fields(model)

struct ConductiveFlux{K}
    conductivity :: K
end

ConductiveFlux(FT::DataType=Float64; conductivity) = ConductiveFlux(convert(FT, conductivity))

@inline function slab_internal_thermal_flux(conductive_flux::ConductiveFlux,
                                            top_temperature,
                                            bottom_temperature,
                                            ice_thickness)

    k = conductive_flux.conductivity
    Tu = top_temperature
    Tb = bottom_temperature
    h = ice_thickness

    return ifelse(h <= 0, zero(h), - k * (Tu - Tb) / h)
end

@inline function slab_internal_thermal_flux(i, j, grid,
                                            top_temperature::Number,
                                            clock, fields, parameters)
    flux = parameters.flux
    bottom_bc = parameters.bottom_thermal_boundary_condition
    liquidus = parameters.liquidus
    Tu = top_temperature
    Tb = bottom_temperature(i, j, grid, bottom_bc, liquidus)
    h = @inbounds fields.h[i, j, 1]
    return slab_internal_thermal_flux(flux, Tu, Tb, h)
end

"""
    SlabSeaIceModel(grid; kw...)

Pretty simple model for sea ice.
"""
function SlabSeaIceModel(grid;
                         clock                             = Clock{eltype(grid)}(0, 0, 1),
                         ice_thickness                     = Field{Center, Center, Nothing}(grid),
                         ice_salinity                      = 0, # psu
                         top_temperature                   = nothing,
                         top_thermal_flux                  = 0,
                         bottom_thermal_flux               = 0,
                         top_thermal_boundary_condition    = MeltingConstrainedFluxBalance(),
                         bottom_thermal_boundary_condition = IceWaterThermalEquilibrium(),
                         internal_thermal_flux             = ConductiveFlux(eltype(grid), conductivity=2),
                         phase_transitions                 = PhaseTransitions(eltype(grid)))

    # Only one time-stepper is supported currently
    timestepper = ForwardEulerTimestepper()
    FT = eltype(grid)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?

    # Wrap ice_salinity in a field
    ice_salinity = field((Center, Center, Nothing), ice_salinity, grid)

    # Construct default top temperature if one is not provided
    if isnothing(top_temperature)
        # Check top boundary condition
        if top_thermal_boundary_condition isa PrescribedTemperature  
            top_temperature = top_thermal_boundary_condition.temperature 
        else # build the default
            top_temperature = Field{Center, Center, Nothing}(grid)
        end
    end

    # Convert to `field` (does nothing if it's already a Field)
    top_temperature = field((Center, Center, Nothing), top_temperature, grid)

    # Construct an internal thermal flux function that captures the liquidus and
    # bottom boundary condition.
    parameters = (flux = internal_thermal_flux,
                  liquidus = phase_transitions.liquidus,
                  bottom_thermal_boundary_condition = bottom_thermal_boundary_condition)

    internal_thermal_flux_function = FluxFunction(slab_internal_thermal_flux;
                                                  parameters,
                                                  top_temperature_dependent=true)

    # Package the external fluxes and boundary conditions
    external_thermal_fluxes = (top = top_thermal_flux,    
                               bottom = bottom_thermal_flux) 

    thermal_boundary_conditions = (top = top_thermal_boundary_condition,
                                   bottom = bottom_thermal_boundary_condition)

    return SlabSeaIceModel(grid,
                           clock,
                           timestepper,
                           ice_thickness,
                           top_temperature,
                           ice_salinity,
                           external_thermal_fluxes,
                           thermal_boundary_conditions,
                           internal_thermal_flux_function,
                           phase_transitions)
end

function set!(model::SSIM; h=nothing)
    !isnothing(h) && set!(model.ice_thickness, h)
    return nothing
end

function time_step!(model::SSIM, Δt; callbacks=nothing)
    grid = model.grid
    Nx, Ny, Nz = size(grid)

    phase_transitions = model.phase_transitions
    liquidus = phase_transitions.liquidus

    h = model.ice_thickness
    Tu = model.top_temperature
    model_fields = fields(model)
    ice_salinity = model.ice_salinity
    clock = model.clock
    bottom_thermal_bc = model.thermal_boundary_conditions.bottom
    top_thermal_bc = model.thermal_boundary_conditions.top

    Qi = model.internal_thermal_flux
    Qb = model.external_thermal_fluxes.bottom
    Qs = model.external_thermal_fluxes.top

    for i = 1:Nx, j=1:Ny

        @inbounds begin
            Tuⁿ = Tu[i, j, 1]
            S = ice_salinity[i, j, 1]
        end

        # 1. Update thickness with ForwardEuler step
        # 1a. Calculate residual fluxes across the top top and bottom
        δQu = top_flux_imbalance(i, j, grid, top_thermal_bc, Tuⁿ, Qi, Qs, clock, model_fields)
        δQb = bottom_flux_imbalance(i, j, grid, bottom_thermal_bc,  Tuⁿ, Qi, Qb, clock, model_fields)

        # 1b. Impose an implicit flux balance if Tu ≤ Tm.
        Tₘ = melting_temperature(liquidus, S) 
        δQu = ifelse(Tuⁿ <= Tₘ, zero(grid), δQu)

        # 1c. Calculate the latent heat of fusion at the top and bottom
        Tbⁿ = bottom_temperature(i, j, grid, bottom_thermal_bc, liquidus)
        ℰu = latent_heat(phase_transitions, Tuⁿ)
        ℰb = latent_heat(phase_transitions, Tbⁿ)

        # 1d. Increment the thickness
        @inbounds begin
            h⁺ = h[i, j, 1] + Δt * (δQb / ℰb + δQu / ℰu)
            h⁺ = max(zero(grid), h⁺)
            h[i, j, 1] = h⁺
        end

        if !isa(top_thermal_bc, PrescribedTemperature)
            # 2. Update top temperature
            Tu⁺ = top_temperature(i, j, grid, top_thermal_bc, Tuⁿ, Qi, Qs, clock, model_fields)

            # 3. Set top temperature (could skip this for `PrescribedTemperature`)
            @inbounds Tu[i, j, 1] = Tu⁺ 
        end
    end

    tick!(model.clock, Δt)

    return nothing
end

function update_state!(model::SSIM)

    top_thermal_bc = model.thermal_boundary_conditions.top

    if !isa(top_thermal_bc, PrescribedTemperature)
        grid = model.grid
        Nx, Ny, Nz = size(grid)

        Tu = model.top_temperature
        model_fields = fields(model)
        clock = model.clock
        top_thermal_bc = model.thermal_boundary_conditions.top
        Qi = model.internal_thermal_flux
        Qs = model.external_thermal_fluxes.top

        # Update top temperature
        for i = 1:Nx, j=1:Ny
            Tu⁻ = @inbounds Tu[i, j, 1]
            Tuⁿ = top_temperature(i, j, grid, top_thermal_bc, Tu⁻, Qi, Qs, clock, model_fields)
            @inbounds Tu[i, j, 1] = Tuⁿ
        end

    end

    return nothing
end

end # module
