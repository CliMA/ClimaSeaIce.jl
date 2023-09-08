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
    top_surface_temperature,
    top_flux_imbalance,
    bottom_flux_imbalance,
    getflux,
    flux_summary

# using RootSolvers: find_zero

using Oceananigans.Architectures: architecture
using Oceananigans.Fields: Field, Center, ZeroField, ConstantField
using Oceananigans.TimeSteppers: Clock, tick!
using Oceananigans.Utils: prettysummary

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

struct SlabSeaIceModel{GR, CL, TS, IT, IC, ST, IS, STF, TBC, CF, P, MIT}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # State
    ice_thickness :: IT
    ice_concentration :: IC
    top_surface_temperature :: ST
    ice_salinity :: IS
    # Boundary conditions
    external_thermal_fluxes :: STF
    thermal_boundary_conditions :: TBC
    # Internal flux
    internal_thermal_flux :: CF
    # Melting and freezing stuff
    phase_transitions :: P
    ice_consolidation_thickness :: MIT
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
    print(io, "├── minimium_ice_thickness: ", prettysummary(model.ice_consolidation_thickness), '\n')
    print(io, "└── external_thermal_fluxes: ", '\n')
    print(io, "    ├── top: ", flux_summary(model.external_thermal_fluxes.top, "    │"), '\n')
    print(io, "    └── bottom: ", flux_summary(model.external_thermal_fluxes.bottom, "     "))
end
         
initialize!(::SSIM) = nothing
default_included_properties(::SSIM) = tuple(:grid)

fields(model::SSIM) = (h = model.ice_thickness,
                       α = model.ice_concentration,
                       Tᵤ = model.top_surface_temperature,
                       Sᵢ = model.ice_salinity)

# TODO: make this correct
prognostic_fields(model::SSIM) = fields(model)

struct ConductiveFlux{K}
    conductivity :: K
end

ConductiveFlux(FT::DataType=Float64; conductivity) = ConductiveFlux(convert(FT, conductivity))

@inline function slab_internal_thermal_flux(conductive_flux::ConductiveFlux,
                                            top_surface_temperature,
                                            bottom_temperature,
                                            ice_thickness)

    k = conductive_flux.conductivity
    Tu = top_surface_temperature
    Tb = bottom_temperature
    h = ice_thickness

    return ifelse(h <= 0, zero(h), - k * (Tu - Tb) / h)
end

@inline function slab_internal_thermal_flux(i, j, grid,
                                            top_surface_temperature::Number,
                                            clock, fields, parameters)
    flux = parameters.flux
    bottom_bc = parameters.bottom_thermal_boundary_condition
    liquidus = parameters.liquidus
    Tu = top_surface_temperature
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
                         ice_consolidation_thickness       = 0.0, # m
                         ice_concentration                 = Field{Center, Center, Nothing}(grid),
                         ice_salinity                      = 0, # psu
                         top_surface_temperature           = nothing,
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
    ice_consolidation_thickness = field((Center, Center, Nothing), ice_consolidation_thickness, grid)

    # Construct default top temperature if one is not provided
    if isnothing(top_surface_temperature)
        # Check top boundary condition
        if top_thermal_boundary_condition isa PrescribedTemperature  
            top_surface_temperature = top_thermal_boundary_condition.temperature 
        else # build the default
            top_surface_temperature = Field{Center, Center, Nothing}(grid)
        end
    end

    # Convert to `field` (does nothing if it's already a Field)
    top_surface_temperature = field((Center, Center, Nothing), top_surface_temperature, grid)

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
                           ice_concentration,
                           top_surface_temperature,
                           ice_salinity,
                           external_thermal_fluxes,
                           thermal_boundary_conditions,
                           internal_thermal_flux_function,
                           phase_transitions,
                           ice_consolidation_thickness)
end

function set!(model::SSIM; h=nothing, α=nothing)
    !isnothing(h) && set!(model.ice_thickness, h)
    !isnothing(α) && set!(model.ice_conentration, α)
    return nothing
end

function time_step!(model::SSIM, Δt; callbacks=nothing)
    grid = model.grid
    Nx, Ny, Nz = size(grid)

    phase_transitions = model.phase_transitions
    liquidus = phase_transitions.liquidus

    h = model.ice_thickness
    hᶜ = model.ice_consolidation_thickness
    α = model.ice_concentration

    Tu = model.top_surface_temperature
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
            hⁿ = h[i, j, 1]
            h★ = hᶜ[i, j, 1]    
            Tuⁿ = Tu[i, j, 1]
            Su = ice_salinity[i, j, 1]
        end

        # Thickness change due to accretion and melting, restricted by minimum allowable value
        consolidated_ice = hⁿ >= h★

        # Compute bottom temperature
        Tbⁿ = bottom_temperature(i, j, grid, bottom_thermal_bc, liquidus)

        if !isa(top_thermal_bc, PrescribedTemperature) # update surface temperature?
            if consolidated_ice # slab is consolidated and has an independent surface temperature
                Tu⁻ = Tuⁿ
                Tuⁿ = top_surface_temperature(i, j, grid, top_thermal_bc, Tu⁻, Qi, Qs, clock, model_fields)
            else # slab is unconsolidated and does not have an independent surface temperature
                Tuⁿ = Tbⁿ
            end

            @inbounds Tu[i, j, 1] = Tuⁿ
        end
        
        # 1. Update thickness with ForwardEuler step
        # 1a. Calculate residual fluxes across the top top and bottom
        δQu = top_flux_imbalance(i, j, grid, top_thermal_bc, Tuⁿ, Qi, Qs, clock, model_fields)

        # δQb = bottom_flux_imbalance(i, j, grid, bottom_thermal_bc, Tuⁿ, Qi, Qb, clock, model_fields)
        #
        # Here we
        #   - distinguish between frazil ice formation (when Qb > 0) and accretion (due to Qi < 0)
        #   - model increases in the ice concentration due to frazil ice formation (rather than simply
        #     modeling frazil ice formation as an increase in _ice thickness_, ie identically as
        #     accretion)
        
        Qiᵢ = getflux(Qi, i, j, grid, Tuⁿ, clock, model_fields)
        Qbᵢ = getflux(Qb, i, j, grid, Tuⁿ, clock, model_fields)

        # 1b. Impose an implicit flux balance if Tu ≤ Tm.
        Tuₘ = melting_temperature(liquidus, Su) # at the top
        δQu = ifelse(Tuⁿ <= Tuₘ, zero(grid), δQu)

        # 1c. Calculate the latent heat of fusion at the top and bottom
        ℰu = latent_heat(phase_transitions, Tuⁿ)
        ℰb = latent_heat(phase_transitions, Tbⁿ)

        # 1d. Increment the thickness
        @inbounds begin
            # Neglect conductive and surface fluxes if ice is below minimum thickness
            Δh = consolidated_ice * Δt * (Qiᵢ / ℰb + δQu / ℰu)

            # "Add" frazil contribution to ice growth (ice formation involves heat flux
            # from ice to ocean thus Qbᵢ < 0):
            Δh -= Δt * Qbᵢ / ℰb

            # Update ice thickness, clipping negative values
            h⁺ = hⁿ + Δh
            h[i, j, 1] = max(zero(grid), h⁺)

            # That's certainly a simple model for ice concentration
            α[i, j, 1] = ifelse(consolidated_ice, one(grid), zero(grid))
        end
    end

    tick!(model.clock, Δt)

    return nothing
end

function update_state!(model::SSIM)

    #=
    top_thermal_bc = model.thermal_boundary_conditions.top

    if !isa(top_thermal_bc, PrescribedTemperature)
        grid = model.grid
        Nx, Ny, Nz = size(grid)

        Tu = model.top_surface_temperature
        model_fields = fields(model)
        clock = model.clock
        top_thermal_bc = model.thermal_boundary_conditions.top
        Qi = model.internal_thermal_flux
        Qs = model.external_thermal_fluxes.top

        # Update top temperature
        for i = 1:Nx, j=1:Ny
            Tu⁻ = @inbounds Tu[i, j, 1]
            Tuⁿ = top_surface_temperature(i, j, grid, top_thermal_bc, Tu⁻, Qi, Qs, clock, model_fields)
            @inbounds Tu[i, j, 1] = Tuⁿ
        end

    end
    =#

    return nothing
end

end # module
