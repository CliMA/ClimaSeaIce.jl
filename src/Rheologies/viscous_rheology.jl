using Base: Callable
using Oceananigans.Fields: AbstractField, location
using Oceananigans.TurbulenceClosures: convert_diffusivity

struct ViscousRheology{N}
    ν :: N
end

const c = Center()
const f = Face()

@inline νᶠᶜᶜ(i, j, k, grid, loc, ν::Number, args...) = ν
@inline νᶜᶠᶠ(i, j, k, grid, loc, ν::Number, args...) = ν

# Array / Field at `Center, Center, Center`
const Lᶜᶜᵃ = Tuple{Center, Center, <:Any}
@inline νᶠᶜᶜ(i, j, k, grid, ::Lᶜᶜᵃ, ν::AbstractArray, args...) = ℑxᵃᶠᵃ(i, j, k, grid, ν)
@inline νᶜᶠᶜ(i, j, k, grid, ::Lᶜᶜᵃ, ν::AbstractArray, args...) = ℑyᵃᶠᵃ(i, j, k, grid, ν)

@inline νᶠᶜᶜ(i, j, k, grid, loc, ν::Callable, clock, args...) = ν(node(i, j, k, grid, f, c, c)..., clock.time)
@inline νᶜᶠᶜ(i, j, k, grid, loc, ν::Callable, clock, args...) = ν(node(i, j, k, grid, c, f, c)..., clock.time)

ViscousRheology(FT::DataType=Float64; ν = 1000.0) = ViscousRheology(convert_diffusivity(FT, ν))

@inline viscosity_location(ν) = (Center(), Center(), Center())
@inline viscosity_location(ν::AbstractField) = location(ν)

@inline ice_stress_ux(i, j, k, grid, vr::ViscousRheology, clock, fields) = νᶜᶠᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * δxᶜᶠᶜ(i, j, k, grid, fields.u)
@inline ice_stress_vx(i, j, k, grid, vr::ViscousRheology, clock, fields) = νᶜᶠᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * δxᶜᶠᶜ(i, j, k, grid, fields.v)
@inline ice_stress_uy(i, j, k, grid, vr::ViscousRheology, clock, fields) = νᶠᶜᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * δyᶠᶜᶜ(i, j, k, grid, fields.u)
@inline ice_stress_vy(i, j, k, grid, vr::ViscousRheology, clock, fields) = νᶠᶜᶜ(i, j, k, grid, viscosity_location(vr.ν), vr.ν, clock, fields) * δyᶠᶜᶜ(i, j, k, grid, fields.v)
