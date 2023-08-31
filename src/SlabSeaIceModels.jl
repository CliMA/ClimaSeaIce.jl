module SlabSeaIceModels

using ClimaSeaIce:
    LinearLiquidus,
    ForwardEulerTimestepper,
    melting_temperature,
    reference_temperature

using ClimaSeaIce.ThermalBoundaryConditions:
    MeltingConstrainedFluxBalance,
    IceWaterThermalEquilibrium,
    PrescribedTemperature,
    FluxFunction,
    SurfaceTemperatureDependent,
    bottom_temperature,
    top_temperature,
    top_flux_imbalance,
    bottom_flux_imbalance

# using RootSolvers: find_zero

using Oceananigans.Fields: Field, Center, ZeroField, ConstantField
using Oceananigans.TimeSteppers: Clock, tick!

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: field, set!
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!, initialize!
import Oceananigans.OutputWriters: default_included_properties

# TODO: move to Oceananigans
field(loc, a::Number, grid) = ConstantField(a)

struct SlabSeaIceModel{GR, CL, TS, IT, ST, IS, STF, TBC, CF, FT, LI}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # State
    ice_thickness :: IT
    surface_temperature :: ST
    ice_salinity :: IS
    # Boundary conditions
    external_thermal_fluxes :: STF
    thermal_boundary_conditions :: TBC
    # Internal flux
    internal_thermal_flux :: CF
    # Ice properties
    ice_density :: FT
    reference_fusion_enthalpy :: FT
    reference_temperature :: FT
    liquidus :: LI
end

const SSIM = SlabSeaIceModel

initialize!(::SSIM) = nothing
default_included_properties(::SSIM) = tuple(:grid)

fields(model::SSIM) = (h = model.ice_thickness,
                       T = model.surface_temperature,
                       S = model.ice_salinity)

# TODO: make this correct
prognostic_fields(model::SSIM) = fields(model)


struct ConductiveFlux{K}
    conductivity :: K
end

ConductiveFlux(FT::DataType=Float64; conductivity) = ConductiveFlux(convert(FT, conductivity))

@inline function slab_internal_thermal_flux(conductive_flux::ConductiveFlux,
                                            surface_temperature,
                                            bottom_temperature,
                                            ice_thickness)

    k = conductive_flux.conductivity
    Tu = surface_temperature
    Tb = bottom_temperature
    h = ice_thickness

    return - k * (Tu - Tb) / h
end

@inline function slab_internal_thermal_flux(i, j, grid,
                                            surface_temperature::Number,
                                            clock, fields, parameters)
    flux = parameters.flux
    bottom_bc = parameters.bottom_thermal_boundary_condition
    liquidus = parameters.liquidus
    Tu = surface_temperature
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
                         surface_temperature               = Field{Center, Center, Nothing}(grid),
                         top_external_thermal_fluxes       = tuple(),
                         bottom_external_thermal_fluxes    = tuple(),
                         top_thermal_boundary_condition    = MeltingConstrainedFluxBalance(),
                         bottom_thermal_boundary_condition = IceWaterThermalEquilibrium(),
                         internal_thermal_flux             = ConductiveFlux(eltype(grid), conductivity=2),
                         ice_density                       = 917,    # kg m⁻³
                         reference_fusion_enthalpy         = 334e3,  # J kg⁻¹
                         reference_temperature             = 273.15, # K
                         liquidus                          = LinearLiquidus())

    timestepper = ForwardEulerTimestepper()
    FT = eltype(grid)

    ice_salinity = field((Center, Center, Nothing), ice_salinity, grid)
    surface_temperature = field((Center, Center, Nothing), surface_temperature, grid)

    parameters = (flux = internal_thermal_flux,
                  liquidus = liquidus,
                  bottom_thermal_boundary_condition = bottom_thermal_boundary_condition)

    internal_thermal_flux_function = FluxFunction(slab_internal_thermal_flux;
                                                  parameters,
                                                  surface_temperature_dependent=true)

    external_thermal_fluxes = (top = top_external_thermal_fluxes,    
                               bottom = bottom_external_thermal_fluxes) 

    thermal_boundary_conditions = (top = top_thermal_boundary_condition,
                                   bottom = bottom_thermal_boundary_condition)

    return SlabSeaIceModel(grid,
                           clock,
                           timestepper,
                           ice_thickness,
                           surface_temperature,
                           ice_salinity,
                           external_thermal_fluxes,
                           thermal_boundary_conditions,
                           internal_thermal_flux_function,
                           convert(FT, ice_density),
                           convert(FT, reference_fusion_enthalpy),
                           convert(FT, reference_temperature),
                           liquidus)
end

# Qs = C * ρₐ * cₐ * u★ * (θₐ - θᵢ)


function time_step!(model::SSIM, Δt; callbacks=nothing)
    grid = model.grid
    Nx, Ny, Nz = size(grid)

    ρᵢ = model.ice_density
    ℒ₀ = model.reference_fusion_enthalpy
    h = model.ice_thickness
    Tu = model.surface_temperature

    model_fields = fields(model)
    ice_salinity = model.ice_salinity
    liquidus = model.liquidus
    clock = model.clock
    bottom_thermal_bc = model.thermal_boundary_conditions.bottom
    top_thermal_bc = model.thermal_boundary_conditions.top

    Qi = model.internal_thermal_flux
    Qb = model.external_thermal_fluxes.bottom
    Qu = model.external_thermal_fluxes.top

    for i = 1:Nx, j=1:Ny

        # 1. Update thickness with ForwardEuler step
        Tuⁿ = @inbounds Tu[i, j, 1]
        δQu =    top_flux_imbalance(i, j, grid, top_thermal_bc,    Tuⁿ, Qi, Qu, clock, model_fields)
        δQb = bottom_flux_imbalance(i, j, grid, bottom_thermal_bc, Tuⁿ, Qi, Qb, clock, model_fields)

        # Impose an implicit flux balance if Tu ≤ Tm.
        S = ice_salinity[i, j, 1]
        Tₘ = melting_temperature(liquidus, S) 
        δQu = ifelse(Tuⁿ <= Tₘ, zero(grid), δQu)

        @inbounds begin
            # TODO: use temperature-specific ℒ rather than reference value
            h⁺ = h[i, j, 1] + Δt * (δQb / (ρᵢ * ℒ₀) + δQu / (ρᵢ * ℒ₀))
            h⁺ = max(zero(grid), h⁺)
            h[i, j, 1] = h⁺
        end

        if !isa(top_thermal_bc, PrescribedTemperature)
            # 2. Update surface temperature
            Tu⁺ = top_temperature(i, j, grid, top_thermal_bc, Tuⁿ, Qi, Qu, clock, model_fields)

            # 3. Set surface temperature (could skip this for `PrescribedTemperature`)
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

        Tu = model.surface_temperature
        model_fields = fields(model)
        clock = model.clock
        top_thermal_bc = model.thermal_boundary_conditions.top
        Qi = model.internal_thermal_flux
        Qu = model.external_thermal_fluxes.top

        # Update surface temperature
        for i = 1:Nx, j=1:Ny
            Tu⁻ = @inbounds Tu[i, j, 1]
            Tuⁿ = top_temperature(i, j, grid, top_thermal_bc, Tu⁻, Qi, Qu, clock, model_fields)
            @inbounds Tu[i, j, 1] = Tuⁿ
        end

    end

    return nothing
end

    end # module
