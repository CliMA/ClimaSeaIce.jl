using Oceananigans.Fields: AbstractField, location
using Oceananigans.TurbulenceClosures: νᶜᶜᶜ, νᶠᶜᶠ, νᶜᶠᶠ, νᶠᶠᶜ, convert_diffusivity

struct ViscousRheology{N}
    ν :: N
end

ViscousRheology(FT::DataType=Float64; ν = 1000.0) = ViscousRheology(convert_diffusivity(FT, ν))

@inline viscosity_location(ν) = (Center(), Center(), Center())
@inline viscosity_location(ν::AbstractField) = location(ν)

@inline viscous_flux_ux(i, j, k, grid, vr::ViscousRheology, clock, fields) = - νᶜᶜᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * ∂xᶜᶜᶜ(i, j, k, grid, fields.u)
@inline viscous_flux_vx(i, j, k, grid, vr::ViscousRheology, clock, fields) = - νᶠᶠᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * ∂xᶠᶠᶜ(i, j, k, grid, fields.v)
@inline viscous_flux_uy(i, j, k, grid, vr::ViscousRheology, clock, fields) = - νᶠᶠᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * ∂yᶠᶠᶜ(i, j, k, grid, fields.u)
@inline viscous_flux_vy(i, j, k, grid, vr::ViscousRheology, clock, fields) = - νᶜᶜᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * ∂yᶜᶜᶜ(i, j, k, grid, fields.v)

@inline function ∂ⱼ_σ₁ⱼ(i, j, k, grid, vr::ViscousRheology, clock, fields)
    return 1 / Azᶠᶜᶜ(i, j, k, grid) * (δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶜᶜ, viscous_flux_ux, vr, clock, fields) +
                                       δyᵃᶜᵃ(i, j, k, grid, Δx_qᶠᶠᶜ, viscous_flux_uy, vr, clock, fields))
end

@inline function ∂ⱼ_σ₂ⱼ(i, j, k, grid, vr::ViscousRheology, clock, fields)
    return 1 / Azᶠᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶠᶜ, viscous_flux_vx, vr, clock, fields) +
                                       δyᵃᶠᵃ(i, j, k, grid, Δx_qᶜᶜᶜ, viscous_flux_vy, vr, clock, fields))
end

