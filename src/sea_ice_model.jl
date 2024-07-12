using Oceananigans.Fields: TracerFields
using ClimaSeaIce.SeaIceThermodynamics: themodynamics_consistent_top_heat_flux

struct SeaIceModel{GR, CL, TS, U, T, IT, IC, TD, D, STF, TBC, SMS, A} <: AbstractModel{TS}
    grid :: GR
    clock :: CL
    timestepper :: TS
    # Prognostic State
    velocities :: U
    tracers :: T
    ice_thickness :: IT
    ice_concentration :: IC
    # Thermodynamics
    sea_ice_thermodynamics :: TD
    # Dynamics
    sea_ice_dynamics :: D
    # External boundary conditions
    external_heat_fluxes :: STF
    external_momentum_stresses :: SMS
    # Numerics
    advection :: A
end

function SeaIceModel(grid;
                     clock                  = Clock{eltype(grid)}(time = 0),
                     ice_thickness          = Field{Center, Center, Nothing}(grid),
                     ice_concentration      = Field{Center, Center, Nothing}(grid),
                     ice_salinity           = 0, # psu
                     top_heat_flux          = nothing,
                     bottom_heat_flux       = 0,
                     velocities             = nothing,
                     advection              = nothing,
                     top_momentum_stress    = nothing, # Fix when introducing dynamics
                     tracers                = (),
                     boundary_conditions    = NamedTuple(),
                     sea_ice_thermodynamics = SlabSeaIceThermodynamics(grid),
                     sea_ice_dynamics       = nothing)

    if isnothing(velocities)
        velocities = (u = ZeroField(), v=ZeroField(), w=ZeroField())
    end

    # Only one time-stepper is supported currently
    timestepper = ForwardEulerTimestepper()

    tracers = tupleit(tracers) # supports tracers=:c keyword argument (for example)
    tracers = TracerFields(tracers, grid, boundary_conditions)

    # TODO: pass `clock` into `field`, so functions can be time-dependent?
    # Wrap ice_salinity in a field
    ice_salinity = field((Center, Center, Nothing), ice_salinity, grid)

    top_heat_flux = thermodynamically_consistent_top_heat_flux(top_heat_flux, sea_ice_thermodynamics)

    # Package the external fluxes and boundary conditions
    external_heat_fluxes = (top = top_heat_flux,    
                            bottom = bottom_heat_flux) 

    return SlabSeaIceThermodynamics(grid,
                                    clock,
                                    timestepper,
                                    velocities,
                                    tracers,
                                    ice_thickness,
                                    ice_concentration,
                                    sea_ice_thermodynamics,
                                    sea_ice_dynamics,
                                    external_heat_fluxes,
                                    top_momentum_stress,
                                    advection)
end

function set!(model::SSIM; h=nothing, ℵ=nothing)
    !isnothing(h) && set!(model.ice_thickness, h)
    !isnothing(ℵ) && set!(model.ice_conentration, ℵ)
    return nothing
end

