using Adapt: Adapt
using Oceananigans.Advection: FluxFormAdvection

"""
    PositivityPreservingAdvection(scheme)

Wrap the flux-form advection `scheme` (for example `WENO()`) with a flux limiter that keeps the ice
content `𝓋 = ℵ h`, the concentration `ℵ` and the snow content `𝓋s = ℵ hs` non-negative without clipping.

High-order reconstructions can ask a nearly empty cell to export more ice than it holds. The step then
clips the result to zero, which creates ice. Here, each step first computes the total area, ice and snow
that leave every cell. It then scales the outgoing fluxes of a cell by

```math
r = \\min\\left(1, \\frac{ℵⁿ}{Δt \\, \\mathrm{Out}_ℵ}, \\frac{𝓋ⁿ}{Δt \\, \\mathrm{Out}_𝓋}, \\frac{𝓋sⁿ}{Δt \\, \\mathrm{Out}_{𝓋s}}\\right),
```

where `ⁿ` is the state at the start of the step (or of the Runge--Kutta substep). Each face flux is
scaled by the `r` of its donor cell, so the fluxes still telescope and the total content is conserved.
The concentration is reconstructed with `scheme`. Each face carries the start-of-step thickness of its
donor cell, so the new thickness of a cell is a weighted average of the thicknesses it keeps and receives.
"""
struct PositivityPreservingAdvection{S, R}
    scheme :: S
    outflow_scaling :: R
end

PositivityPreservingAdvection(scheme) = PositivityPreservingAdvection(scheme, nothing)

Base.summary(a::PositivityPreservingAdvection) = string("PositivityPreservingAdvection(", summary(a.scheme), ")")
Base.show(io::IO, a::PositivityPreservingAdvection) = print(io, summary(a))

Oceananigans.Grids.required_halo_size_x(a::PositivityPreservingAdvection) = Oceananigans.Grids.required_halo_size_x(a.scheme)
Oceananigans.Grids.required_halo_size_y(a::PositivityPreservingAdvection) = Oceananigans.Grids.required_halo_size_y(a.scheme)

Adapt.adapt_structure(to, a::PositivityPreservingAdvection) =
    PositivityPreservingAdvection(Adapt.adapt(to, a.scheme), Adapt.adapt(to, a.outflow_scaling))

function Oceananigans.Advection.materialize_advection(a::PositivityPreservingAdvection{<:Any, Nothing}, grid)
    scheme = Oceananigans.Advection.materialize_advection(a.scheme, grid)
    return PositivityPreservingAdvection(scheme, Field{Center, Center, Nothing}(grid))
end

# The limiter keeps the content non-negative, so only round-off concentrations are emptied.
@inline minimum_ice_concentration(FT, ::PositivityPreservingAdvection) = convert(FT, 1e-16)

@inline x_advection(scheme) = scheme
@inline y_advection(scheme) = scheme
@inline x_advection(scheme::FluxFormAdvection) = scheme.x
@inline y_advection(scheme::FluxFormAdvection) = scheme.y

# The state the step (or substep) updates from.
start_of_step_state(model, ::ForwardEulerTimeStepper) = (model.ice_thickness, model.ice_concentration, model.snow_thickness)

function start_of_step_state(model, timestepper)
    Ψ⁻ = timestepper.Ψ⁻
    hsⁿ = isnothing(model.snow_thickness) ? nothing : Ψ⁻.hs
    return (Ψ⁻.h, Ψ⁻.ℵ, hsⁿ)
end

function compute_tracer_tendencies!(model::SIM, advection::PositivityPreservingAdvection, Δt)
    grid = model.grid
    arch = architecture(grid)

    scheme = advection.scheme
    r = advection.outflow_scaling

    hⁿ, ℵⁿ, hsⁿ = start_of_step_state(model, model.timestepper)

    launch!(arch, grid, :xy, _compute_outflow_scaling!, r, grid, scheme, model.velocities,
            model.ice_concentration, ℵⁿ, hⁿ, hsⁿ, Δt)

    # A face flux is scaled by its donor cell, which lies in the halo at the domain edge.
    fill_halo_regions!(r)

    launch!(arch, grid, :xy, _compute_limited_tracer_tendencies!, model.timestepper.Gⁿ, grid, scheme, r,
            model.velocities, model.ice_concentration, hⁿ, hsⁿ)

    return nothing
end

@inline function cell_outflow(i, j, k, grid, flux_x, flux_y, scheme, U, args...)
    Fʷ = flux_2D_x(i,   j, k, grid, flux_x, x_advection(scheme), U.u, args...)
    Fᵉ = flux_2D_x(i+1, j, k, grid, flux_x, x_advection(scheme), U.u, args...)
    Fˢ = flux_2D_y(i, j,   k, grid, flux_y, y_advection(scheme), U.v, args...)
    Fᴺ = flux_2D_y(i, j+1, k, grid, flux_y, y_advection(scheme), U.v, args...)
    outflow = max(Fᵉ, zero(Fᵉ)) + max(-Fʷ, zero(Fʷ)) + max(Fᴺ, zero(Fᴺ)) + max(-Fˢ, zero(Fˢ))
    return outflow / Azᶜᶜᶜ(i, j, k, grid)
end

# Fraction of the outgoing amount `Δc` that the available amount `c` can supply.
@inline available_fraction(c, Δc) = ifelse(Δc > c, max(c, zero(c)) / Δc, one(c))

# The donor cell of a face is set by the sign of the area flux, which can oppose the velocity where WENO undershoots.
@inline function upwind_thickness_flux_x(i, j, k, grid, scheme, u, ℵ, h)
    F = _area_flux_x(i, j, k, grid, scheme, u, ℵ)
    return @inbounds F * ifelse(F > 0, h[i-1, j, k], h[i, j, k])
end

@inline function upwind_thickness_flux_y(i, j, k, grid, scheme, v, ℵ, h)
    F = _area_flux_y(i, j, k, grid, scheme, v, ℵ)
    return @inbounds F * ifelse(F > 0, h[i, j-1, k], h[i, j, k])
end

@inline snow_fraction(i, j, k, grid, scheme, U, ℵ, ℵⁿ, ::Nothing, Δt) = one(grid)

@inline function snow_fraction(i, j, k, grid, scheme, U, ℵ, ℵⁿ, hsⁿ, Δt)
    𝓋sⁿ = @inbounds ℵⁿ[i, j, 1] * hsⁿ[i, j, 1]
    Δ𝓋s = Δt * cell_outflow(i, j, k, grid, upwind_thickness_flux_x, upwind_thickness_flux_y, scheme, U, ℵ, hsⁿ)
    return available_fraction(𝓋sⁿ, Δ𝓋s)
end

@kernel function _compute_outflow_scaling!(r, grid, scheme, U, ℵ, ℵⁿ, hⁿ, hsⁿ, Δt)
    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)

    @inbounds begin
        𝓋ⁿ = ℵⁿ[i, j, 1] * hⁿ[i, j, 1]
        Δ𝓋 = Δt * cell_outflow(i, j, kᴺ, grid, upwind_thickness_flux_x, upwind_thickness_flux_y, scheme, U, ℵ, hⁿ)
        Δℵ = Δt * cell_outflow(i, j, kᴺ, grid, _area_flux_x, _area_flux_y, scheme, U, ℵ)

        rᵢ = min(available_fraction(𝓋ⁿ, Δ𝓋), available_fraction(ℵⁿ[i, j, 1], Δℵ))
        r[i, j, 1] = min(rᵢ, snow_fraction(i, j, kᴺ, grid, scheme, U, ℵ, ℵⁿ, hsⁿ, Δt))
    end
end

@inline function limited_flux_x(i, j, k, grid, flux, scheme, u, r, args...)
    F = flux(i, j, k, grid, scheme, u, args...)
    return @inbounds F * ifelse(F > 0, r[i-1, j, 1], r[i, j, 1])
end

@inline function limited_flux_y(i, j, k, grid, flux, scheme, v, r, args...)
    F = flux(i, j, k, grid, scheme, v, args...)
    return @inbounds F * ifelse(F > 0, r[i, j-1, 1], r[i, j, 1])
end

@inline function limited_divergence(i, j, k, grid, flux_x, flux_y, scheme, U, r, args...)
    return 1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, flux_2D_x, limited_flux_x, flux_x, x_advection(scheme), U.u, r, args...) +
                                       δyᵃᶜᵃ(i, j, k, grid, flux_2D_y, limited_flux_y, flux_y, y_advection(scheme), U.v, r, args...))
end

@kernel function _compute_limited_tracer_tendencies!(Gⁿ, grid, scheme, r, U, ℵ, hⁿ, hsⁿ)
    i, j = @index(Global, NTuple)
    kᴺ = size(grid, 3)

    @inbounds begin
        Gⁿ.h[i, j, 1] = - limited_divergence(i, j, kᴺ, grid, upwind_thickness_flux_x, upwind_thickness_flux_y, scheme, U, r, ℵ, hⁿ)
        Gⁿ.ℵ[i, j, 1] = - limited_divergence(i, j, kᴺ, grid, _area_flux_x, _area_flux_y, scheme, U, r, ℵ)
    end

    compute_limited_snow_tendency!(i, j, kᴺ, Gⁿ, grid, scheme, r, U, ℵ, hsⁿ)
end

@inline compute_limited_snow_tendency!(i, j, k, Gⁿ, grid, scheme, r, U, ℵ, ::Nothing) = nothing

@inline function compute_limited_snow_tendency!(i, j, k, Gⁿ, grid, scheme, r, U, ℵ, hsⁿ)
    @inbounds Gⁿ.hs[i, j, 1] = - limited_divergence(i, j, k, grid, upwind_thickness_flux_x, upwind_thickness_flux_y, scheme, U, r, ℵ, hsⁿ)
    return nothing
end
