using Oceananigans.Advection: FluxFormAdvection,
                              WENO,
                              FunctionStencil,
                              AbstractUpwindBiasedAdvectionScheme,
                              AbstractCenteredAdvectionScheme,
                              bias,
                              _advective_tracer_flux_x,
                              _advective_tracer_flux_y,
                              _biased_interpolate_xᶠᵃᵃ,
                              _biased_interpolate_yᵃᶠᵃ,
                              _symmetric_interpolate_xᶠᵃᵃ,
                              _symmetric_interpolate_yᵃᶠᵃ,
                              conditional_flux_fcc,
                              conditional_flux_cfc

@inline ice_content(i, j, k, grid, ℵ, h) = @inbounds ℵ[i, j, k] * h[i, j, k]
@inline limiter(ϕ, a, b) = clamp(ϕ, min(a, b), max(a, b))

@inline function content_and_area_x(i, j, k, grid, scheme::WENO, b, ℵ, h)
    ss = FunctionStencil(ice_content)
    𝓋ᶠ = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, b, ice_content, ss, ℵ, h)
    ℵᶠ = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, b, ℵ, ss, ℵ, h)
    return 𝓋ᶠ, ℵᶠ
end

@inline function content_and_area_y(i, j, k, grid, scheme::WENO, b, ℵ, h)
    ss = FunctionStencil(ice_content)
    𝓋ᶠ = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, b, ice_content, ss, ℵ, h)
    ℵᶠ = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, b, ℵ, ss, ℵ, h)
    return 𝓋ᶠ, ℵᶠ
end

@inline function content_and_area_x(i, j, k, grid, scheme::AbstractCenteredAdvectionScheme, b, ℵ, h)
    𝓋ᶠ = _symmetric_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, b, ice_content, ℵ, h)
    ℵᶠ = _symmetric_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, b, ℵ)
    return 𝓋ᶠ, ℵᶠ
end

@inline function content_and_area_y(i, j, k, grid, scheme::AbstractCenteredAdvectionScheme, b, ℵ, h)
    𝓋ᶠ = _symmetric_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, b, ice_content, ℵ, h)
    ℵᶠ = _symmetric_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, b, ℵ)
    return 𝓋ᶠ, ℵᶠ
end

@inline function content_and_area_x(i, j, k, grid, scheme::AbstractUpwindBiasedAdvectionScheme, b, ℵ, h)
    𝓋ᶠ = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, b, ice_content, ℵ, h)
    ℵᶠ = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, b, ℵ)
    return 𝓋ᶠ, ℵᶠ
end

@inline function content_and_area_y(i, j, k, grid, scheme::AbstractUpwindBiasedAdvectionScheme, b, ℵ, h)
    𝓋ᶠ = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, b, ice_content, ℵ, h)
    ℵᶠ = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, b, ℵ)
    return 𝓋ᶠ, ℵᶠ
end

@inline function reconstruct_thickness_x(i, j, k, grid, scheme, U, ℵ, h)
    @inbounds begin
        b = bias(U[i, j, k])
        𝓋ᶠ, ℵᶠ = content_and_area_x(i, j, k, grid, scheme, b, ℵ, h)
        hᶠ = ifelse(ℵᶠ > 0, 𝓋ᶠ / ℵᶠ, zero(grid))
        hᶠ = limiter(hᶠ, h[i-1, j, k], h[i, j, k])
    end
    return hᶠ
end

@inline function reconstruct_thickness_y(i, j, k, grid, scheme, V, ℵ, h)
    @inbounds begin
        b = bias(V[i, j, k])
        𝓋ᶠ, ℵᶠ = content_and_area_y(i, j, k, grid, scheme, b, ℵ, h)
        hᶠ = ifelse(ℵᶠ > 0, 𝓋ᶠ / ℵᶠ, zero(grid))
        hᶠ = limiter(hᶠ, h[i, j-1, k], h[i, j, k])
    end
    return hᶠ
end

@inline _advective_thickness_flux_x(i, j, k, grid, scheme, U, ℵ, h) = advective_thickness_flux_x(i, j, k, grid, scheme, U, ℵ, h)
@inline _advective_thickness_flux_y(i, j, k, grid, scheme, U, ℵ, h) = advective_thickness_flux_y(i, j, k, grid, scheme, U, ℵ, h)

@inline _advective_thickness_flux_x(i, j, k, ibg::ImmersedBoundaryGrid, scheme, U, ℵ, h) =
    conditional_flux_fcc(i, j, k, ibg, zero(ibg), advective_thickness_flux_x(i, j, k, ibg, scheme, U, ℵ, h))

@inline _advective_thickness_flux_y(i, j, k, ibg::ImmersedBoundaryGrid, scheme, U, ℵ, h) =
    conditional_flux_cfc(i, j, k, ibg, zero(ibg), advective_thickness_flux_y(i, j, k, ibg, scheme, U, ℵ, h))

@inline reconstruct_area_x(i, j, k, grid, scheme::AbstractUpwindBiasedAdvectionScheme, U, ℵ) = _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, bias(@inbounds U[i, j, k]), ℵ)
@inline reconstruct_area_y(i, j, k, grid, scheme::AbstractUpwindBiasedAdvectionScheme, V, ℵ) = _biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, bias(@inbounds V[i, j, k]), ℵ)
@inline reconstruct_area_x(i, j, k, grid, scheme::AbstractCenteredAdvectionScheme, U, ℵ) = _symmetric_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, ℵ)
@inline reconstruct_area_y(i, j, k, grid, scheme::AbstractCenteredAdvectionScheme, V, ℵ) = _symmetric_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, ℵ)

@inline function area_flux_x(i, j, k, grid, scheme, U, ℵ)
    ℵᴿ = reconstruct_area_x(i, j, k, grid, scheme, U, ℵ)
    @inbounds ℵᶠ = limiter(ℵᴿ, ℵ[i-1, j, k], ℵ[i, j, k])
    return @inbounds Axᶠᶜᶜ(i, j, k, grid) * U[i, j, k] * ℵᶠ
end

@inline function area_flux_y(i, j, k, grid, scheme, V, ℵ)
    ℵᴿ = reconstruct_area_y(i, j, k, grid, scheme, V, ℵ)
    @inbounds ℵᶠ = limiter(ℵᴿ, ℵ[i, j-1, k], ℵ[i, j, k])
    return @inbounds Ayᶜᶠᶜ(i, j, k, grid) * V[i, j, k] * ℵᶠ
end

@inline _area_flux_x(i, j, k, grid, scheme, U, ℵ) = area_flux_x(i, j, k, grid, scheme, U, ℵ)
@inline _area_flux_y(i, j, k, grid, scheme, V, ℵ) = area_flux_y(i, j, k, grid, scheme, V, ℵ)

@inline _area_flux_x(i, j, k, ibg::ImmersedBoundaryGrid, scheme, U, ℵ) = conditional_flux_fcc(i, j, k, ibg, zero(ibg), area_flux_x(i, j, k, ibg, scheme, U, ℵ))
@inline _area_flux_y(i, j, k, ibg::ImmersedBoundaryGrid, scheme, V, ℵ) = conditional_flux_cfc(i, j, k, ibg, zero(ibg), area_flux_y(i, j, k, ibg, scheme, V, ℵ))

# Content flux `Ax·u·ℵ_face · h_face`: the monotone area flux times the limited thickness reconstruction.
@inline advective_thickness_flux_x(i, j, k, grid, advection, U, ℵ, h) =
    area_flux_x(i, j, k, grid, advection, U, ℵ) * reconstruct_thickness_x(i, j, k, grid, advection, U, ℵ, h)

@inline advective_thickness_flux_y(i, j, k, grid, advection, V, ℵ, h) =
    area_flux_y(i, j, k, grid, advection, V, ℵ) * reconstruct_thickness_y(i, j, k, grid, advection, V, ℵ, h)

# Divide by Δzᶠᶜᶜ in case of a moving grid (zstar for example)
@inline function _advective_thickness_flux_2D_x(i, j, k, grid, scheme, U, ℵ, h)
    Δz = Δzᶠᶜᶜ(i, j, k, grid)
    return ifelse(Δz > 0, _advective_thickness_flux_x(i, j, k, grid, scheme, U, ℵ, h) / Δz, zero(grid))
end

@inline function _advective_thickness_flux_2D_y(i, j, k, grid, scheme, U, ℵ, h)
    Δz = Δzᶜᶠᶜ(i, j, k, grid)
    return ifelse(Δz > 0, _advective_thickness_flux_y(i, j, k, grid, scheme, U, ℵ, h) / Δz, zero(grid))
end

@inline function _area_flux_2D_x(i, j, k, grid, scheme, U, ℵ)
    Δz = Δzᶠᶜᶜ(i, j, k, grid)
    return ifelse(Δz > 0, _area_flux_x(i, j, k, grid, scheme, U, ℵ) / Δz, zero(grid))
end

@inline function _area_flux_2D_y(i, j, k, grid, scheme, V, ℵ)
    Δz = Δzᶜᶠᶜ(i, j, k, grid)
    return ifelse(Δz > 0, _area_flux_y(i, j, k, grid, scheme, V, ℵ) / Δz, zero(grid))
end

@inline function _advective_tracer_flux_2D_x(i, j, k, grid, scheme, U, c)
    Δz = Δzᶠᶜᶜ(i, j, k, grid)
    return ifelse(Δz > 0, _advective_tracer_flux_x(i, j, k, grid, scheme, U, c) / Δz, zero(grid))
end

@inline function _advective_tracer_flux_2D_y(i, j, k, grid, scheme, U, c)
    Δz = Δzᶜᶠᶜ(i, j, k, grid)
    return ifelse(Δz > 0, _advective_tracer_flux_y(i, j, k, grid, scheme, U, c) / Δz, zero(grid))
end

@inline div_Uℵh(i, j, k, grid, ::Nothing, U, ℵ, h) = zero(grid)

# Volume-per-area tendency: returns ∇·(U·ℵ·h) with a 2D (area-based) metric.
@inline function div_Uℵh(i, j, k, grid, advection, U, ℵ, h)
    return 1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_thickness_flux_2D_x, advection, U.u, ℵ, h) +
                                       δyᵃᶜᵃ(i, j, k, grid, _advective_thickness_flux_2D_y, advection, U.v, ℵ, h))
end

@inline function div_Uℵh(i, j, k, grid, advection::FluxFormAdvection, U, ℵ, h)
    return 1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_thickness_flux_2D_x, advection.x, U.u, ℵ, h) +
                                       δyᵃᶜᵃ(i, j, k, grid, _advective_thickness_flux_2D_y, advection.y, U.v, ℵ, h))
end

# Concentration tendency: ∇·(U·ℵ) with the monotone area flux and a 2D (area-based) metric.
@inline div_Uℵ(i, j, k, grid, ::Nothing, U, ℵ) = zero(grid)

@inline function div_Uℵ(i, j, k, grid, advection, U, ℵ)
    return 1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _area_flux_2D_x, advection, U.u, ℵ) +
                                       δyᵃᶜᵃ(i, j, k, grid, _area_flux_2D_y, advection, U.v, ℵ))
end

@inline function div_Uℵ(i, j, k, grid, advection::FluxFormAdvection, U, ℵ)
    return 1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _area_flux_2D_x, advection.x, U.u, ℵ) +
                                       δyᵃᶜᵃ(i, j, k, grid, _area_flux_2D_y, advection.y, U.v, ℵ))
end
