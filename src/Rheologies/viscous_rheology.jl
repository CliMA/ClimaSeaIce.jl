using Oceananigans.Fields: AbstractField, location
using Oceananigans.TurbulenceClosures: νᶜᶜᶜ, νᶠᶜᶠ, νᶜᶠᶠ, νᶠᶠᶜ, convert_diffusivity

struct ViscousRheology{N}
    ν :: N
end

ViscousRheology(FT::DataType=Float64; ν = 1000.0) = ViscousRheology(convert_diffusivity(FT, ν))

@inline viscosity_location(ν) = (Center(), Center(), Center())
@inline viscosity_location(ν::AbstractField) = location(ν)

@inline ice_stress_ux(i, j, k, grid, vr::ViscousRheology, clock, fields) = νᶜᶠᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * δxᶜᶠᶜ(i, j, k, grid, fields.u)
@inline ice_stress_vx(i, j, k, grid, vr::ViscousRheology, clock, fields) = νᶜᶠᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * δxᶜᶠᶜ(i, j, k, grid, fields.v)
@inline ice_stress_uy(i, j, k, grid, vr::ViscousRheology, clock, fields) = νᶜᶠᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * δyᶠᶜᶜ(i, j, k, grid, fields.u)
@inline ice_stress_vy(i, j, k, grid, vr::ViscousRheology, clock, fields) = νᶜᶠᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * δyᶠᶜᶜ(i, j, k, grid, fields.v)
