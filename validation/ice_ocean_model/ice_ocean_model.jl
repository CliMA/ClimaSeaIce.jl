using Oceananigans.Operators

using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Models: AbstractModel
using Oceananigans.TimeSteppers: tick!
using Oceananigans.Utils: launch!

using KernelAbstractions: @kernel, @index
using KernelAbstractions.Extras.LoopInfo: @unroll

using ClimaOcean
using ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: compute_atmosphere_ocean_fluxes,
                                                     compute_atmosphere_sea_ice_fluxes,
                                                     compute_sea_ice_ocean_salinity_flux!,
                                                     sea_ice_ocean_latent_heat_flux!

import Oceananigans.TimeSteppers: update_state!
import ClimaOcean.OceanSeaIceModels.CrossRealmFluxes: compute_sea_ice_ocean_fluxes

function update_state!(coupled_model::OceanSeaIceModel, callbacks=[]; compute_tendencies=false)
    time = Time(coupled_model.clock.time)
    update_model_field_time_series!(coupled_model.atmosphere, time)
    compute_atmosphere_ocean_fluxes!(coupled_model) 
    compute_sea_ice_ocean_fluxes!(coupled_model)
    # compute_atmosphere_sea_ice_fluxes!(coupled_model)
    return nothing
end

function compute_sea_ice_ocean_fluxes!(coupled_model)
    compute_sea_ice_ocean_salinity_flux!(coupled_model)
    sea_ice_ocean_latent_heat_flux!(coupled_model)
    return nothing
end
