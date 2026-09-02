using Adapt: Adapt
using Oceananigans.Grids: halo_size
using Oceananigans.ImmersedBoundaries: immersed_cell
using Oceananigans.Operators: Δxᶜᶜᶜ, Δyᶜᶜᶜ, Δxᶠᶜᶜ, Δyᶠᶜᶜ, Δxᶜᶠᶜ, Δyᶜᶠᶜ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ
using Oceananigans.Utils: KernelParameters

"""
    struct IncrementalRemapping

Two-dimensional conservative transport of the ice concentration `ℵ` and of the thickness-like tracers
`h` and `hs`, following [Lipscomb and Hunke (2004)](@cite LipscombHunke2004).

`ℵ`, `h` and `hs` are reconstructed as limited linear functions inside each cell. The material crossing
a face during a step occupies the region swept by the backward trajectories of the face, and both the
transported area `𝔉_ℵ = ∫_R ℵ̃ dA` and the transported content `𝔉_𝓋 = ∫_R ℵ̃ h̃ dA` are integrals of the
same reconstruction over that same region, which keeps `h = 𝓋/ℵ` inside the range of `h̃` at any `ℵ`.
The integrals are evaluated by tensor Gauss--Legendre quadrature over the bilinear parametrization of
the swept region.

The transport integrates over the whole step, so the scheme requires `timestepper = :ForwardEuler`, a
Courant number below one in each direction, and a halo of at least two.

Keyword arguments
=================

- `quadrature_nodes`: number of Gauss--Legendre nodes per parametric direction. Default: `3`.
"""
struct IncrementalRemapping{N, R, T}
    reconstruction :: R
    transports :: T
end

IncrementalRemapping(; quadrature_nodes::Int = 3) = IncrementalRemapping{quadrature_nodes, Nothing, Nothing}(nothing, nothing)

Base.summary(::IncrementalRemapping{N}) where N = string("IncrementalRemapping(quadrature_nodes=", N, ")")
Base.show(io::IO, scheme::IncrementalRemapping) = print(io, summary(scheme))

Oceananigans.Grids.required_halo_size_x(::IncrementalRemapping) = 2
Oceananigans.Grids.required_halo_size_y(::IncrementalRemapping) = 2

function Adapt.adapt_structure(to, scheme::IncrementalRemapping{N}) where N
    reconstruction = Adapt.adapt(to, scheme.reconstruction)
    transports = Adapt.adapt(to, scheme.transports)
    return IncrementalRemapping{N, typeof(reconstruction), typeof(transports)}(reconstruction, transports)
end

function Oceananigans.Advection.materialize_advection(scheme::IncrementalRemapping{N, Nothing}, grid) where N
    Hx, Hy, _ = halo_size(grid)

    if Hx < 2 || Hy < 2
        throw(ArgumentError("IncrementalRemapping requires a horizontal halo of at least 2, got ($Hx, $Hy)"))
    end

    center_field() = Field{Center, Center, Nothing}(grid)

    # `(xa, ya)` is the centroid of the reconstructed ice area, which anchors the tracer reconstructions.
    reconstruction = (ℵx = center_field(), ℵy = center_field(),
                      hx = center_field(), hy = center_field(),
                      hsx = center_field(), hsy = center_field(),
                      xa = center_field(), ya = center_field())

    transports = (x = (ℵ = Field{Face, Center, Nothing}(grid),
                       𝓋 = Field{Face, Center, Nothing}(grid),
                       𝓋s = Field{Face, Center, Nothing}(grid)),
                  y = (ℵ = Field{Center, Face, Nothing}(grid),
                       𝓋 = Field{Center, Face, Nothing}(grid),
                       𝓋s = Field{Center, Face, Nothing}(grid)))

    R = typeof(reconstruction)
    T = typeof(transports)

    return IncrementalRemapping{N, R, T}(reconstruction, transports)
end

validate_advection_timestepper(::IncrementalRemapping, timestepper) =
    throw(ArgumentError("IncrementalRemapping requires timestepper = :ForwardEuler, got $(summary(timestepper))"))

validate_advection_timestepper(::IncrementalRemapping, ::ForwardEulerTimeStepper) = nothing

@inline minimum_ice_concentration(FT, ::IncrementalRemapping) = zero(FT)

#####
##### Gauss--Legendre rules on the unit interval
#####

@inline quadrature_rule(FT, ::Val{1}) = ((FT(1/2), FT(1)),)

@inline quadrature_rule(FT, ::Val{2}) = ((FT(1/2) - FT(0.28867513459481287), FT(1/2)),
                                         (FT(1/2) + FT(0.28867513459481287), FT(1/2)))

@inline quadrature_rule(FT, ::Val{3}) = ((FT(1/2) - FT(0.3872983346207417), FT(5/18)),
                                         (FT(1/2),                          FT(8/18)),
                                         (FT(1/2) + FT(0.3872983346207417), FT(5/18)))

@inline quadrature_rule(FT, ::Val{4}) = ((FT(1/2) - FT(0.43056815579702629), FT(0.17392742256872693)),
                                         (FT(1/2) - FT(0.16999052179242816), FT(0.32607257743127307)),
                                         (FT(1/2) + FT(0.16999052179242816), FT(0.32607257743127307)),
                                         (FT(1/2) + FT(0.43056815579702629), FT(0.17392742256872693)))

#####
##### Limited linear reconstruction
#####

@inline submerged(i, j, k, grid) = false
@inline submerged(i, j, k, ibg::ImmersedBoundaryGrid) = immersed_cell(i, j, k, ibg)

@inline wet_value(i, j, k, grid, c, c₀) = @inbounds ifelse(submerged(i, j, k, grid), c₀, c[i, j, k])

@inline function neighborhood_extrema(i, j, k, grid, c)
    @inbounds c₀ = c[i, j, k]

    cᵂ  = wet_value(i-1, j,   k, grid, c, c₀)
    cᴱ  = wet_value(i+1, j,   k, grid, c, c₀)
    cˢ  = wet_value(i,   j-1, k, grid, c, c₀)
    cᴺ  = wet_value(i,   j+1, k, grid, c, c₀)
    cˢᵂ = wet_value(i-1, j-1, k, grid, c, c₀)
    cˢᴱ = wet_value(i+1, j-1, k, grid, c, c₀)
    cᴺᵂ = wet_value(i-1, j+1, k, grid, c, c₀)
    cᴺᴱ = wet_value(i+1, j+1, k, grid, c, c₀)

    cmin = min(c₀, cᵂ, cᴱ, cˢ, cᴺ, cˢᵂ, cˢᴱ, cᴺᵂ, cᴺᴱ)
    cmax = max(c₀, cᵂ, cᴱ, cˢ, cᴺ, cˢᵂ, cˢᴱ, cᴺᴱ, cᴺᵂ)
    return cmin, cmax
end

# Thickness of a neighbour, replaced by the local value where that neighbour carries no ice, so open
# water does not drag the thickness reconstruction of an ice-covered cell toward zero.
@inline function iced_thickness(i, j, k, grid, h, ℵ, h₀)
    @inbounds ℵᵢ = ℵ[i, j, k]
    @inbounds hᵢ = h[i, j, k]
    ℵᵐⁱⁿ = minimum_ice_concentration(typeof(h₀))
    iced = (ℵᵢ > ℵᵐⁱⁿ) & !submerged(i, j, k, grid)
    return ifelse(iced, hᵢ, h₀)
end

@inline function iced_extrema(i, j, k, grid, h, ℵ, h₀)
    hᵂ  = iced_thickness(i-1, j,   k, grid, h, ℵ, h₀)
    hᴱ  = iced_thickness(i+1, j,   k, grid, h, ℵ, h₀)
    hˢ  = iced_thickness(i,   j-1, k, grid, h, ℵ, h₀)
    hᴺ  = iced_thickness(i,   j+1, k, grid, h, ℵ, h₀)
    hˢᵂ = iced_thickness(i-1, j-1, k, grid, h, ℵ, h₀)
    hˢᴱ = iced_thickness(i+1, j-1, k, grid, h, ℵ, h₀)
    hᴺᵂ = iced_thickness(i-1, j+1, k, grid, h, ℵ, h₀)
    hᴺᴱ = iced_thickness(i+1, j+1, k, grid, h, ℵ, h₀)

    hmin = min(h₀, hᵂ, hᴱ, hˢ, hᴺ, hˢᵂ, hˢᴱ, hᴺᵂ, hᴺᴱ)
    hmax = max(h₀, hᵂ, hᴱ, hˢ, hᴺ, hˢᵂ, hˢᴱ, hᴺᵂ, hᴺᴱ)
    return hmin, hmax
end

# Barth--Jespersen scaling: the largest α ≤ 1 keeping the reconstruction within the neighbourhood
# bounds over a cell whose extreme corner sits `δ` away from the anchor.
@inline function limiting_factor(δ, c₀, cmin, cmax)
    α = min(one(δ), (cmax - c₀) / δ, (c₀ - cmin) / δ)
    return ifelse(δ > 0, α, one(δ))
end

@inline function centered_x_gradient(i, j, k, grid, c)
    @inbounds c₀ = c[i, j, k]
    δc = wet_value(i+1, j, k, grid, c, c₀) - wet_value(i-1, j, k, grid, c, c₀)
    return δc / (Δxᶠᶜᶜ(i, j, k, grid) + Δxᶠᶜᶜ(i+1, j, k, grid))
end

@inline function centered_y_gradient(i, j, k, grid, c)
    @inbounds c₀ = c[i, j, k]
    δc = wet_value(i, j+1, k, grid, c, c₀) - wet_value(i, j-1, k, grid, c, c₀)
    return δc / (Δyᶜᶠᶜ(i, j, k, grid) + Δyᶜᶠᶜ(i, j+1, k, grid))
end

@inline function iced_x_gradient(i, j, k, grid, h, ℵ, h₀)
    δh = iced_thickness(i+1, j, k, grid, h, ℵ, h₀) - iced_thickness(i-1, j, k, grid, h, ℵ, h₀)
    return δh / (Δxᶠᶜᶜ(i, j, k, grid) + Δxᶠᶜᶜ(i+1, j, k, grid))
end

@inline function iced_y_gradient(i, j, k, grid, h, ℵ, h₀)
    δh = iced_thickness(i, j+1, k, grid, h, ℵ, h₀) - iced_thickness(i, j-1, k, grid, h, ℵ, h₀)
    return δh / (Δyᶜᶠᶜ(i, j, k, grid) + Δyᶜᶠᶜ(i, j+1, k, grid))
end

# Anchoring at the area centroid `(xa, ya)` makes the area-weighted mean of the reconstruction return
# the cell mean exactly.
@inline function limited_tracer_gradients(i, j, k, grid, h, ℵ, xa, ya, Δx, Δy)
    @inbounds h₀ = h[i, j, k]

    hx = iced_x_gradient(i, j, k, grid, h, ℵ, h₀)
    hy = iced_y_gradient(i, j, k, grid, h, ℵ, h₀)

    hmin, hmax = iced_extrema(i, j, k, grid, h, ℵ, h₀)

    δh = abs(hx) * (Δx / 2 + abs(xa)) + abs(hy) * (Δy / 2 + abs(ya))
    α = limiting_factor(δh, h₀, hmin, hmax)

    return α * hx, α * hy
end

@inline limited_tracer_gradients(i, j, k, grid, ::Nothing, ℵ, xa, ya, Δx, Δy) = zero(grid), zero(grid)

@kernel function _compute_remapping_reconstruction!(reconstruction, grid, ℵ, h, hs)
    i, j = @index(Global, NTuple)
    k = size(grid, 3)

    Δx = Δxᶜᶜᶜ(i, j, k, grid)
    Δy = Δyᶜᶜᶜ(i, j, k, grid)

    @inbounds ℵ₀ = ℵ[i, j, k]

    ℵx = centered_x_gradient(i, j, k, grid, ℵ)
    ℵy = centered_y_gradient(i, j, k, grid, ℵ)

    ℵmin, ℵmax = neighborhood_extrema(i, j, k, grid, ℵ)

    δℵ = abs(ℵx) * Δx / 2 + abs(ℵy) * Δy / 2
    α = limiting_factor(δℵ, ℵ₀, ℵmin, ℵmax)

    α = ifelse(submerged(i, j, k, grid), zero(α), α)

    ℵx = α * ℵx
    ℵy = α * ℵy

    # xa = ∫ℵ̃ x dA / ∫ℵ̃ dA = ℵx Δx² / (12 ℵ₀), measured from the cell centre, and likewise in y
    ℵᵐⁱⁿ = minimum_ice_concentration(typeof(ℵ₀))
    iced = ℵ₀ > ℵᵐⁱⁿ
    xa = ifelse(iced, ℵx * Δx^2 / (12 * ℵ₀), zero(ℵ₀))
    ya = ifelse(iced, ℵy * Δy^2 / (12 * ℵ₀), zero(ℵ₀))

    hx, hy = limited_tracer_gradients(i, j, k, grid, h, ℵ, xa, ya, Δx, Δy)
    hsx, hsy = limited_tracer_gradients(i, j, k, grid, hs, ℵ, xa, ya, Δx, Δy)

    @inbounds begin
        reconstruction.ℵx[i, j, 1] = ℵx
        reconstruction.ℵy[i, j, 1] = ℵy
        reconstruction.hx[i, j, 1] = hx
        reconstruction.hy[i, j, 1] = hy
        reconstruction.hsx[i, j, 1] = hsx
        reconstruction.hsy[i, j, 1] = hsy
        reconstruction.xa[i, j, 1] = xa
        reconstruction.ya[i, j, 1] = ya
    end
end

#####
##### Point evaluation of the reconstruction
#####

# Evaluate the reconstructions at `(x, y)`, measured from the centre of cell `(i, j)`. The point may lie
# in a neighbouring cell, so the containing cell is located first and its own reconstruction is used.
@inline function reconstruct(i, j, k, grid, x, y, ℵ, h, hs, r)
    Δx = Δxᶜᶜᶜ(i, j, k, grid)
    Δy = Δyᶜᶜᶜ(i, j, k, grid)

    di = ifelse(x < -Δx/2, -1, ifelse(x > Δx/2, 1, 0))
    dj = ifelse(y < -Δy/2, -1, ifelse(y > Δy/2, 1, 0))

    I = i + di
    J = j + dj

    xc = ifelse(di == -1, -Δxᶠᶜᶜ(i, j, k, grid), ifelse(di == 1, Δxᶠᶜᶜ(i+1, j, k, grid), zero(grid)))
    yc = ifelse(dj == -1, -Δyᶜᶠᶜ(i, j, k, grid), ifelse(dj == 1, Δyᶜᶠᶜ(i, j+1, k, grid), zero(grid)))

    δx = x - xc
    δy = y - yc

    @inbounds begin
        ℵ̃ = ℵ[I, J, k] + r.ℵx[I, J, 1] * δx + r.ℵy[I, J, 1] * δy
        δxa = δx - r.xa[I, J, 1]
        δya = δy - r.ya[I, J, 1]
    end

    h̃ = reconstruct_tracer(I, J, k, h, r.hx, r.hy, δxa, δya)
    h̃s = reconstruct_tracer(I, J, k, hs, r.hsx, r.hsy, δxa, δya)

    return max(ℵ̃, zero(ℵ̃)), h̃, h̃s
end

@inline reconstruct_tracer(I, J, k, ::Nothing, hx, hy, δx, δy) = zero(δx)
@inline reconstruct_tracer(I, J, k, h, hx, hy, δx, δy) = @inbounds h[I, J, k] + hx[I, J, 1] * δx + hy[I, J, 1] * δy

#####
##### Transport through the swept region of a face
#####

# The region swept through the face is spanned by (x, y)(η, s) = p(η) - s Δt 𝐮(η) on (η, s) ∈ [0, 1]²,
# where `p(η)` runs along the face and `𝐮(η)` is the corner velocity interpolated along it, in
# coordinates local to the face midpoint. The Jacobian is signed, so a face whose normal velocity
# changes sign along its length contributes the difference of the two lobes.
@inline function x_face_transport(i, j, k, grid, ::Val{N}, Δt, u, v, ℵ, h, hs, r) where N
    FT = eltype(grid)

    Δy = Δyᶠᶜᶜ(i, j, k, grid)
    xᶜ = Δxᶜᶜᶜ(i, j, k, grid) / 2

    u₁ = ℑyᵃᶠᵃ(i, j,   k, grid, u)
    v₁ = ℑxᶠᵃᵃ(i, j,   k, grid, v)
    u₂ = ℑyᵃᶠᵃ(i, j+1, k, grid, u)
    v₂ = ℑxᶠᵃᵃ(i, j+1, k, grid, v)

    δu = u₂ - u₁
    δv = v₂ - v₁

    𝔉ℵ = zero(FT)
    𝔉𝓋 = zero(FT)
    𝔉𝓈 = zero(FT)

    rule = quadrature_rule(FT, Val(N))

    for (η, wη) in rule, (s, ws) in rule
        uη = u₁ + η * δu
        vη = v₁ + η * δv

        x = - s * Δt * uη
        y = - Δy/2 + η * Δy - s * Δt * vη

        𝒥 = Δt * uη * Δy + s * Δt^2 * (δu * vη - uη * δv)

        ℵ̃, h̃, h̃s = reconstruct(i, j, k, grid, x - xᶜ, y, ℵ, h, hs, r)

        w = wη * ws * 𝒥
        𝔉ℵ += w * ℵ̃
        𝔉𝓋 += w * ℵ̃ * h̃
        𝔉𝓈 += w * ℵ̃ * h̃s
    end

    return 𝔉ℵ, 𝔉𝓋, 𝔉𝓈
end

@inline function y_face_transport(i, j, k, grid, ::Val{N}, Δt, u, v, ℵ, h, hs, r) where N
    FT = eltype(grid)

    Δx = Δxᶜᶠᶜ(i, j, k, grid)
    yᶜ = Δyᶜᶜᶜ(i, j, k, grid) / 2

    u₁ = ℑyᵃᶠᵃ(i,   j, k, grid, u)
    v₁ = ℑxᶠᵃᵃ(i,   j, k, grid, v)
    u₂ = ℑyᵃᶠᵃ(i+1, j, k, grid, u)
    v₂ = ℑxᶠᵃᵃ(i+1, j, k, grid, v)

    δu = u₂ - u₁
    δv = v₂ - v₁

    𝔉ℵ = zero(FT)
    𝔉𝓋 = zero(FT)
    𝔉𝓈 = zero(FT)

    rule = quadrature_rule(FT, Val(N))

    for (ξ, wξ) in rule, (s, ws) in rule
        uξ = u₁ + ξ * δu
        vξ = v₁ + ξ * δv

        x = - Δx/2 + ξ * Δx - s * Δt * uξ
        y = - s * Δt * vξ

        𝒥 = Δt * vξ * Δx + s * Δt^2 * (δv * uξ - vξ * δu)

        ℵ̃, h̃, h̃s = reconstruct(i, j, k, grid, x, y - yᶜ, ℵ, h, hs, r)

        w = wξ * ws * 𝒥
        𝔉ℵ += w * ℵ̃
        𝔉𝓋 += w * ℵ̃ * h̃
        𝔉𝓈 += w * ℵ̃ * h̃s
    end

    return 𝔉ℵ, 𝔉𝓋, 𝔉𝓈
end

# Nothing crosses an immersed face.
@inline mask_immersed_transport_x(i, j, k, grid, 𝔉) = 𝔉
@inline mask_immersed_transport_y(i, j, k, grid, 𝔉) = 𝔉

@inline mask_immersed_transport_x(i, j, k, ibg::ImmersedBoundaryGrid, 𝔉) =
    conditional_flux_fcc(i, j, k, ibg, zero(ibg), 𝔉)

@inline mask_immersed_transport_y(i, j, k, ibg::ImmersedBoundaryGrid, 𝔉) =
    conditional_flux_cfc(i, j, k, ibg, zero(ibg), 𝔉)

@kernel function _compute_x_transports!(transports, grid, nodes, Δt, u, v, ℵ, h, hs, reconstruction)
    i, j = @index(Global, NTuple)
    k = size(grid, 3)

    𝔉ℵ, 𝔉𝓋, 𝔉𝓈 = x_face_transport(i, j, k, grid, nodes, Δt, u, v, ℵ, h, hs, reconstruction)

    𝔉ℵ = mask_immersed_transport_x(i, j, k, grid, 𝔉ℵ)
    𝔉𝓋 = mask_immersed_transport_x(i, j, k, grid, 𝔉𝓋)
    𝔉𝓈 = mask_immersed_transport_x(i, j, k, grid, 𝔉𝓈)

    @inbounds begin
        transports.ℵ[i, j, 1] = 𝔉ℵ
        transports.𝓋[i, j, 1] = 𝔉𝓋
        transports.𝓋s[i, j, 1] = 𝔉𝓈
    end
end

@kernel function _compute_y_transports!(transports, grid, nodes, Δt, u, v, ℵ, h, hs, reconstruction)
    i, j = @index(Global, NTuple)
    k = size(grid, 3)

    𝔉ℵ, 𝔉𝓋, 𝔉𝓈 = y_face_transport(i, j, k, grid, nodes, Δt, u, v, ℵ, h, hs, reconstruction)

    𝔉ℵ = mask_immersed_transport_y(i, j, k, grid, 𝔉ℵ)
    𝔉𝓋 = mask_immersed_transport_y(i, j, k, grid, 𝔉𝓋)
    𝔉𝓈 = mask_immersed_transport_y(i, j, k, grid, 𝔉𝓈)

    @inbounds begin
        transports.ℵ[i, j, 1] = 𝔉ℵ
        transports.𝓋[i, j, 1] = 𝔉𝓋
        transports.𝓋s[i, j, 1] = 𝔉𝓈
    end
end

#####
##### Tendencies
#####

@inline function transport_divergence(i, j, 𝔉x, 𝔉y)
    @inbounds return 𝔉x[i+1, j, 1] - 𝔉x[i, j, 1] + 𝔉y[i, j+1, 1] - 𝔉y[i, j, 1]
end

@kernel function _compute_remapping_tendencies!(Gⁿ, grid, transports, Δt, snow_thickness)
    i, j = @index(Global, NTuple)
    k = size(grid, 3)

    V = Azᶜᶜᶜ(i, j, k, grid) * Δt

    @inbounds begin
        Gⁿ.ℵ[i, j, 1] = - transport_divergence(i, j, transports.x.ℵ, transports.y.ℵ) / V
        Gⁿ.h[i, j, 1] = - transport_divergence(i, j, transports.x.𝓋, transports.y.𝓋) / V
    end

    compute_remapping_snow_tendency!(i, j, Gⁿ, transports, V, snow_thickness)
end

@inline compute_remapping_snow_tendency!(i, j, Gⁿ, transports, V, ::Nothing) = nothing

@inline function compute_remapping_snow_tendency!(i, j, Gⁿ, transports, V, snow_thickness)
    @inbounds Gⁿ.hs[i, j, 1] = - transport_divergence(i, j, transports.x.𝓋s, transports.y.𝓋s) / V
    return nothing
end

function compute_tracer_tendencies!(model::SIM, advection::IncrementalRemapping{N}, Δt) where N
    grid = model.grid
    arch = architecture(grid)

    Nx, Ny, _ = size(grid)

    ℵ = model.ice_concentration
    h = model.ice_thickness
    hs = model.snow_thickness

    u = model.velocities.u
    v = model.velocities.v

    reconstruction = advection.reconstruction
    transports = advection.transports

    # The transport through a boundary face reads the reconstruction of the cell outside it, so the
    # first halo cell is reconstructed too.
    reconstruction_parameters = KernelParameters(0:Nx+1, 0:Ny+1)
    launch!(arch, grid, reconstruction_parameters, _compute_remapping_reconstruction!, reconstruction, grid, ℵ, h, hs)

    nodes = Val(N)

    launch!(arch, grid, KernelParameters(1:Nx+1, 1:Ny), _compute_x_transports!, transports.x, grid, nodes, Δt, u, v, ℵ, h, hs, reconstruction)
    launch!(arch, grid, KernelParameters(1:Nx, 1:Ny+1), _compute_y_transports!, transports.y, grid, nodes, Δt, u, v, ℵ, h, hs, reconstruction)
    launch!(arch, grid, :xy, _compute_remapping_tendencies!, model.timestepper.Gⁿ, grid, transports, Δt, hs)

    return nothing
end
