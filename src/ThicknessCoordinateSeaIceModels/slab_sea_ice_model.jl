using ClimaSeaIce:
    LinearLiquidus,
    ConstantBulkSalinity,
    ForwardEulerTimestepper,
    melting_temperature,
    reference_temperature

using ClimaSeaIce.IceSurfaces:
    MeltingRadiatingSurface,
    IceWaterSurface,
    TimeLinearizedSurfaceTemperatureSolver

using Oceananigans.Fields: ZeroField

# Simulations interface
import Oceananigans: fields, prognostic_fields
import Oceananigans.Fields: set!
import Oceananigans.TimeSteppers: time_step!, update_state!
import Oceananigans.Simulations: reset!, initialize!

struct SlabSeaIceModel{G, C, TS, F, IS, A, O, FT, LI}
    grid :: G
    clock :: C
    timestepper :: TS
    # Ice state
    ice_thickness :: F
    top_ice_temperature :: F
    ice_salinity :: IS
    # Boundary conditions
    top_surface :: A
    bottom_surface :: O
    # Ice properties
    ice_conductivity :: FT
    ice_density :: FT
    reference_fusion_enthalpy :: FT
    liquidus :: LI
end

const SSIM = SlabSeaIceModel

initialize!(model::SSIM) = nothing

fields(model::SSIM) = (h = model.ice_thickness,
                       T = model.top_ice_temperature,
                       S = model.ice_salinity)

# TODO: make this correct
prognostic_fields(model::SSIM) = fields(model)

"""
    SlabSeaIceModel(grid; kw...)

Pretty simple model for sea ice.
"""
function SlabSeaIceModel(grid;
                         clock = Clock{eltype(grid)}(0, 0, 1),
                         ice_thickness = Field{Center, Center, Nothing}(grid),
                         ice_salinity = ConstantBulkSalinity(eltype(grid)),
                         top_ice_temperature = Field{Center, Center, Nothing}(grid),
                         top_surface = MeltingRadiatingSurface(eltype(grid)),
                         bottom_surface = IceWaterSurface(ZeroField()),
                         ice_conductivity = 2.0, # W m² / ᵒC
                         ice_density = 911.0, # kg m⁻³
                         reference_fusion_enthalpy = 333.55, # kJ kg⁻¹
                         liquidus = LinearLiquidus())

    timestepper = ForwardEulerTimestepper()
    FT = eltype(grid)

    return SlabSeaIceModel(grid,
                           clock,
                           timestepper,
                           ice_thickness,
                           top_ice_temperature,
                           ice_salinity,
                           top_surface,
                           bottom_surface,
                           convert(FT, ice_conductivity),
                           convert(FT, ice_density),
                           convert(FT, reference_fusion_enthalpy),
                           liquidus)
end

const PTS = PrescribedTemperatureSurface
const CBS = ConstantBulkSalinity
const TLMRS = MeltingRadiatingSurface{<:TimeLinearizedSurfaceTemperatureSolver}

@inline top_ice_salinity(i, j, grid, salinity::CBS, fields) = salinity.constant
@inline bottom_surface_temperature(i, j, grid, bs::PTS, liquidus) = bs.temperature[i, j, 1]
@inline bottom_surface_flux_imbalance(::PTS, Qex, k, h, Tu, Tb) = zero(Qex)

@inline function bottom_surface_temperature(i, j, grid, bottom_surface::IceWaterSurface, liquidus)
    Sₒ = @inbounds bottom_surface.ocean_salinity[i, j, 1]
    return melting_temperature(liquidus, Sₒ)
end

@inline bottom_surface_flux_imbalance(::IceWaterSurface, Qex, k, h, Tu, Tb) = Qex - k * (Tu - Tb) / h

@inline function top_surface_flux_imbalance(surface::MeltingRadiatingSurface, Qx, k, h, Tu, Tb)
    σ = surface.stefan_boltzmann_constant
    ϵ = surface.emissivity
    T₀ = surface.reference_temperature

    # Emitted radiation
    Qe = ϵ * σ * (Tu + T₀)^4

    # Conductive flux
    Qk = - k * (Tu - Tb) / h

    # Return the sum of all three
    return Qx + Qe + Qk
end

@inline function top_surface_temperature(surface::TLMRS, Qex, k, h, Tu, Tb)
    δQ = top_surface_flux_imbalance(surface, Qex, k, h, Tu, Tb)

    σ = surface.stefan_boltzmann_constant
    ϵ = surface.emissivity
    T₀ = surface.reference_temperature

    # Derivative of ΣQ with respect to T
    dQ_dT = 4 * ϵ * σ * (Tu + T₀)^3 - k / h

    @show dQ_dT

    ΔT = - δQ / dQ_dT

    return Tu + ΔT
end

function time_step!(model::SSIM, Δt; callbacks=nothing)
    grid = model.grid
    Nx, Ny, Nz = size(grid)

    ρᵢ = model.ice_density
    ℒ₀ = model.reference_fusion_enthalpy
    k = model.ice_conductivity
    h = model.ice_thickness
    Tu = model.top_ice_temperature
    bottom_surface = model.bottom_surface

    top_surface = model.top_surface
    σ = model.top_surface.stefan_boltzmann_constant

    model_fields = fields(model)
    ice_salinity = model.ice_salinity
    liquidus = model.liquidus

    for i = 1:Nx, j=1:Ny
        @inbounds begin 
            Tuᵢ = Tu[i, j, 1]
            hᵢ = h[i, j, 1]
            Qex = top_surface.extrinsic_heat_flux[i, j, 1]
        end

        # Compute bottom temperature (depends on ocean state)
        Tbᵢ = bottom_surface_temperature(i, j, grid, bottom_surface, liquidus, model_fields)

        # Compute surface temperature (depends on bottom temperature and extrinsic heat fluxes)
        Tuᵢ = top_surface_temperature(top_surface, Qex, k, hᵢ, Tuᵢ, Tbᵢ)

        # Clip surface temperature to melting temperature
        Suᵢ = top_ice_salinity(i, j, grid, ice_salinity, model_fields)
        Tmᵢ = melting_temperature(liquidus, Suᵢ)
        Tuᵢ = min(Tmᵢ, Tuᵢ)

        # Recompute surface flux imbalance
        δQu = top_surface_flux_imbalance(top_surface, Qex, k, hᵢ, Tuᵢ, Tbᵢ)
        δQu = ifelse(Tuᵢ == Tmᵢ, δQu, zero(δQu)) # clip
        Δhu = δQu / (ρᵢ * ℒ₀) # approximation for large Stefan number 

        δQb = bottom_surface_flux_imbalance(bottom_surface, Qex, k, hᵢ, Tuᵢ, Tbᵢ)
        Δhb = δQb / (ρᵢ * ℒ₀) # approximation for large Stefan number 

        @inbounds begin
            h[i, j, 1] += (Δhu + Δhb) * Δt
            Tu[i, j, 1] = Tuᵢ
        end
    end

    tick!(model.clock, Δt)

    return nothing
end

update_state!(model::SSIM) = nothing

