using Oceananigans

struct CenteredWENO5 end

  reconstruct_on_node(i, j, k, grid, scheme, σ, args...) = ℑxyᶠᶠᵃ(i, j, k, grid, σ, args...)
reconstruct_on_center(i, j, k, grid, scheme, σ, args...) = ℑxyᶜᶜᵃ(i, j, k, grid, σ, args...)

  reconstruct_on_node(i, j, k, grid, scheme::CenteredWENO5, σ, args...) = interpolate_xyᶠᶠ(i, j, k, grid, scheme, σ, args...)
reconstruct_on_center(i, j, k, grid, scheme::CenteredWENO5, σ, args...) = interpolate_xyᶜᶜ(i, j, k, grid, scheme, σ, args...)

using Oceananigans.Advection: stencil_coefficients, uniform_reconstruction_coefficients

########### Centered WENO Interpolation scheme ############
#
# We use 4 3-point stencils to compute the 5th order WENO interpolation
# This interpolation should limit to a 6th order centered scheme
#
# Computation of the optimal coefficients:
# C6 = uniform_reconstruction_coefficients(Float64, Val(:symmetric), 3)

# Left biased stencils are the S0 and S1
const S0 = ( 2.0, -7.0, 11.0) ./ 6 # stencil_coefficients(Float64, 50, 2, collect(1:100), collect(1:100); order=3)
const S1 = (-1.0,  5.0,  2.0) ./ 6 # stencil_coefficients(Float64, 50, 1, collect(1:100), collect(1:100); order=3)

# Right biased stencils are the S2 and S3
const S2 = reverse(S1)
const S3 = reverse(S0)

# Generally speaking then:
const C★3 = 0.05 # C6[6] / S3[3]
const C★0 = 0.05 # C6[1] / S0[1]

# And the other two coefficients:
const C★1 = 0.45 # (C6[2] - C★0 * S0[2]) / S1[1] 
const C★2 = 0.45 # (C6[5] - C★3 * S3[2]) / S2[3]

@inline   left_biased_β(FT, ψ) = @inbounds FT(13/12) * (ψ[1] - 2ψ[2] + ψ[3])^2 + FT(1/4) * ( ψ[1] - 4ψ[2] + 3ψ[3])^2
@inline center_biased_β(FT, ψ) = @inbounds FT(13/12) * (ψ[1] - 2ψ[2] + ψ[3])^2 + FT(1/4) * ( ψ[1]         -  ψ[3])^2
@inline  right_biased_β(FT, ψ) = @inbounds FT(13/12) * (ψ[1] - 2ψ[2] + ψ[3])^2 + FT(1/4) * (3ψ[1] - 4ψ[2] +  ψ[3])^2

@inline getvalue(i, j, k, grid, σ::Function, args...)      = σ(i, j, k, grid, args...)
@inline getvalue(i, j, k, grid, σ::AbstractArray, args...) = @inbounds σ[i, j, k]

const ϵ = 1e-10

# @inline function centered_weno(σ₀::FT, σ₁, σ₂, σ₃, σ₄, σ₅) where FT
    
#     S₀ = (σ₀, σ₁, σ₂)
#     S₁ = (σ₁, σ₂, σ₃)
#     S₂ = (σ₂, σ₃, σ₄)
#     S₃ = (σ₃, σ₄, σ₅)

#     β₀ =    left_biased_β(FT, S₀)
#     β₁ = (center_biased_β(FT, S₁) +  left_biased_β(FT, S₁)) / 2
#     β₂ = (center_biased_β(FT, S₂) + right_biased_β(FT, S₂)) / 2
#     β₃ =   right_biased_β(FT, S₃)

#     ω₀ = C★0 / (β₀ + ϵ)^2
#     ω₁ = C★1 / (β₁ + ϵ)^2
#     ω₂ = C★2 / (β₂ + ϵ)^2
#     ω₃ = C★3 / (β₃ + ϵ)^2

#     q₀ = sum(S₀ .* S0)
#     q₁ = sum(S₁ .* S1)
#     q₂ = sum(S₂ .* S2)
#     q₃ = sum(S₃ .* S3)
    
#     return (ω₀ * q₀ + ω₁ * q₁ + ω₂ * q₂ + ω₃ * q₃) / (ω₀ + ω₁ + ω₂ + ω₃)
# end

const S02 = ( 1, -5, 13, 3) ./ 12
const S12 = (-1,  5,  2)    ./ 6
const S22 = ( 2,  5, -1)    ./ 6
const S32 = ( 3, 13, -5, 1) ./ 12

const S13 = ( 1, -2, 0,  2, -1) ./ 24
const S23 = (-1,  2, 0, -2,  1) ./ 24

@inline βᵂ₀(S) = (2S[1] - 9S[2] + 18S[3] - 11S[4])^2 / 36 + (1S[1] - 4S[2] + 5S[3] - 2S[4])^2 / 12 * 13 + (1S[1] - 3S[2] + 3S[3] - 1S[4])^2 / 720 * 781
@inline βᵂ₃(S) = (11S[1] - 18S[2] + 9S[3] - 2S[4])^2 / 36 + (2S[1] - 5S[2] + 4S[3] - 1S[4])^2 / 12 * 13 + (- 1S[1] + 3S[2] - 3S[3] + 1S[4])^2 / 720 * 781
@inline βᵂ₁(S) = (1S[1] - 2S[2] + 1S[3])^2 / 12 * 13 + (S[1] - S[3])^2 / 4
@inline βᵂ₂(S) = (1S[1] - 2S[2] + 1S[3])^2 / 12 * 13 + (3S[1] - 4S[2] + S[3])^2 / 4
@inline βᶜ₁(S) = (- 1S[1] + 2S[2] - 2S[4] + 1S[5]) * ((- 1S[1] + 2S[2] - 2S[4] + 1S[5]) / 320 * 89  + (S[2] - S[4]) / 12)
@inline βᶜ₂(S) = (- 1S[1] + 2S[2] - 2S[4] + 1S[5]) * ((- 1S[1] + 2S[2] - 2S[4] + 1S[5]) / 547 * 960 - (19S[3] - 34S[4] + 15S[5]) / 12)

@inline function centered_weno(σ₀::FT, σ₁, σ₂, σ₃, σ₄, σ₅) where FT
    
    S₀ = (σ₀, σ₁, σ₂, σ₃)
    S₁ = (σ₁, σ₂, σ₃)
    S₂ = (σ₂, σ₃, σ₄)
    S₃ = (σ₂, σ₃, σ₄, σ₅)
    Sᶜ = (σ₀, σ₁, σ₂, σ₃, σ₄)

    β₀ = βᵂ₀(S₀) 
    β₁ = βᵂ₁(S₁) + βᶜ₁(Sᶜ) 
    β₂ = βᵂ₂(S₂) + βᶜ₂(Sᶜ) 
    β₃ = βᵂ₃(S₃)  

    ω₀ = 2 / (β₀ + ϵ)^2  # C★0
    ω₁ = 3 / (β₁ + ϵ)^2  # C★1
    ω₂ = 3 / (β₂ + ϵ)^2  # C★2
    ω₃ = 2 / (β₃ + ϵ)^2  # C★3

    q₀ = sum(S₀ .* S02)
    q₁ = sum(S₁ .* S12) + sum(S13 .* Sᶜ)
    q₂ = sum(S₂ .* S22) + sum(S23 .* Sᶜ)
    q₃ = sum(S₃ .* S32)
    
    return (ω₀ * q₀ + ω₁ * q₁ + ω₂ * q₂ + ω₃ * q₃) / (ω₀ + ω₁ + ω₂ + ω₃)
end

# Now we can compute the interpolation:
@inline function interpolate_xᶠ(i, j, k, grid, ::CenteredWENO5, σ, args...)
    σ₀ = getvalue(i-3, j, k, grid, σ, args...)
    σ₁ = getvalue(i-2, j, k, grid, σ, args...)
    σ₂ = getvalue(i-1, j, k, grid, σ, args...)
    σ₃ = getvalue(i  , j, k, grid, σ, args...)
    σ₄ = getvalue(i+1, j, k, grid, σ, args...)
    σ₅ = getvalue(i+2, j, k, grid, σ, args...)
    return centered_weno(σ₀, σ₁, σ₂, σ₃, σ₄, σ₅)
end

# Now we can compute the interpolation:
@inline function interpolate_yᶠ(i, j, k, grid, ::CenteredWENO5, σ, args...)
    σ₀ = getvalue(i, j-3, k, grid, σ, args...)
    σ₁ = getvalue(i, j-2, k, grid, σ, args...)
    σ₂ = getvalue(i, j-1, k, grid, σ, args...)
    σ₃ = getvalue(i, j  , k, grid, σ, args...)
    σ₄ = getvalue(i, j+1, k, grid, σ, args...)
    σ₅ = getvalue(i, j+2, k, grid, σ, args...)
    return centered_weno(σ₀, σ₁, σ₂, σ₃, σ₄, σ₅)
end

# Now we can compute the interpolation:
@inline function interpolate_xᶜ(i, j, k, grid, ::CenteredWENO5, σ, args...)
    σ₀ = getvalue(i-2, j, k, grid, σ, args...)
    σ₁ = getvalue(i-1, j, k, grid, σ, args...)
    σ₂ = getvalue(i  , j, k, grid, σ, args...)
    σ₃ = getvalue(i+1, j, k, grid, σ, args...)
    σ₄ = getvalue(i+2, j, k, grid, σ, args...)
    σ₅ = getvalue(i+3, j, k, grid, σ, args...)
    return centered_weno(σ₀, σ₁, σ₂, σ₃, σ₄, σ₅)
end

# Now we can compute the interpolation:
@inline function interpolate_yᶜ(i, j, k, grid, ::CenteredWENO5, σ, args...)
    σ₀ = getvalue(i, j-2, k, grid, σ, args...)
    σ₁ = getvalue(i, j-1, k, grid, σ, args...)
    σ₂ = getvalue(i, j  , k, grid, σ, args...)
    σ₃ = getvalue(i, j+1, k, grid, σ, args...)
    σ₄ = getvalue(i, j+2, k, grid, σ, args...)
    σ₅ = getvalue(i, j+3, k, grid, σ, args...)
    return centered_weno(σ₀, σ₁, σ₂, σ₃, σ₄, σ₅)
end

function interpolate_xyᶜᶜ(i, j, k, grid, scheme::CenteredWENO5, σ, args...)
    σ₀ = interpolate_xᶜ(i, j-2, k, grid, scheme, σ, args...)
    σ₁ = interpolate_xᶜ(i, j-1, k, grid, scheme, σ, args...)
    σ₂ = interpolate_xᶜ(i, j,   k, grid, scheme, σ, args...)
    σ₃ = interpolate_xᶜ(i, j+1, k, grid, scheme, σ, args...)
    σ₄ = interpolate_xᶜ(i, j+2, k, grid, scheme, σ, args...)
    σ₅ = interpolate_xᶜ(i, j+3, k, grid, scheme, σ, args...)
    return centered_weno(σ₀, σ₁, σ₂, σ₃, σ₄, σ₅)
end

function interpolate_xyᶠᶠ(i, j, k, grid, scheme::CenteredWENO5, σ, args...)
    σ₀ = interpolate_xᶠ(i, j-3, k, grid, scheme, σ, args...)
    σ₁ = interpolate_xᶠ(i, j-2, k, grid, scheme, σ, args...)
    σ₂ = interpolate_xᶠ(i, j-1, k, grid, scheme, σ, args...)
    σ₃ = interpolate_xᶠ(i, j  , k, grid, scheme, σ, args...)
    σ₄ = interpolate_xᶠ(i, j+1, k, grid, scheme, σ, args...)
    σ₅ = interpolate_xᶠ(i, j+2, k, grid, scheme, σ, args...)
    return centered_weno(σ₀, σ₁, σ₂, σ₃, σ₄, σ₅)
end

function interpolate_xyᶠᶜ(i, j, k, grid, scheme::CenteredWENO5, σ, args...)
    σ₀ = interpolate_xᶠ(i, j-2, k, grid, scheme, σ, args...)
    σ₁ = interpolate_xᶠ(i, j-1, k, grid, scheme, σ, args...)
    σ₂ = interpolate_xᶠ(i, j  , k, grid, scheme, σ, args...)
    σ₃ = interpolate_xᶠ(i, j+1, k, grid, scheme, σ, args...)
    σ₄ = interpolate_xᶠ(i, j+2, k, grid, scheme, σ, args...)
    σ₅ = interpolate_xᶠ(i, j+3, k, grid, scheme, σ, args...)
    return centered_weno(σ₀, σ₁, σ₂, σ₃, σ₄, σ₅)
end

function interpolate_xyᶜᶠ(i, j, k, grid, scheme::CenteredWENO5, σ, args...)
    σ₀ = interpolate_xᶜ(i, j-3, k, grid, scheme, σ, args...)
    σ₁ = interpolate_xᶜ(i, j-2, k, grid, scheme, σ, args...)
    σ₂ = interpolate_xᶜ(i, j-1, k, grid, scheme, σ, args...)
    σ₃ = interpolate_xᶜ(i, j  , k, grid, scheme, σ, args...)
    σ₄ = interpolate_xᶜ(i, j+1, k, grid, scheme, σ, args...)
    σ₅ = interpolate_xᶜ(i, j+2, k, grid, scheme, σ, args...)
    return centered_weno(σ₀, σ₁, σ₂, σ₃, σ₄, σ₅)
end

@inline function strain_rate_xy_centered(i, j, k, grid, scheme::CenteredWENO5, u, v) 
    δv = δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, interpolate_xyᶠᶜ, scheme, v) 
    δu = δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, interpolate_xyᶜᶠ, scheme, u)
    return (δv + δu) / Azᶜᶜᶜ(i, j, k, grid) / 2
end