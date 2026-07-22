using ClimaSeaIce
using ClimaSeaIce.Rheologies: ∂ⱼ_σ₁ⱼ, ∂ⱼ_σ₂ⱼ,
                              strain_rate_xx, strain_rate_yy, strain_rate_xy,
                              ElastoViscoPlasticRheology
using Oceananigans
using Oceananigans.Operators
using Oceananigans.Grids: Center, Face, λnodes, φnodes
using Oceananigans.Fields: location
using Test

#####
##### Discrete energy budget for the internal-stress operators.
#####
##### The volume-integrated work done by the stress divergence on the velocity field
##### must equal minus the volume-integrated stress power σ:ε̇ (up to boundary fluxes):
#####
#####   ∑ [u ∂ⱼσ₁ⱼ + v ∂ⱼσ₂ⱼ] Az  =  - ∑ [σ₁₁ ε̇₁₁ + σ₂₂ ε̇₂₂ + 2 σ₁₂ ε̇₁₂] Az  +  (boundary)
#####
##### This is a summation-by-parts identity: it holds to machine precision iff `∂ⱼσ` is the
##### negative adjoint of the discrete strain-rate operator under the area-weighted inner
##### product. The metric-aware ("Δ²-inside") formulation is an adjoint pair with the
##### strain rate and closes the budget on a curvilinear grid; the plain flux-form
##### divergence used previously is not the adjoint there and leaks energy.

# Previous (flux-form) stress divergence, reproduced here for contrast.
@inline Δyᶜᶜ_σ(i, j, k, grid, σ) = @inbounds Δyᶜᶜᶜ(i, j, k, grid) * σ[i, j, k]
@inline Δxᶠᶠ_σ(i, j, k, grid, σ) = @inbounds Δxᶠᶠᶜ(i, j, k, grid) * σ[i, j, k]
@inline Δyᶠᶠ_σ(i, j, k, grid, σ) = @inbounds Δyᶠᶠᶜ(i, j, k, grid) * σ[i, j, k]
@inline Δxᶜᶜ_σ(i, j, k, grid, σ) = @inbounds Δxᶜᶜᶜ(i, j, k, grid) * σ[i, j, k]

@inline old_∂ⱼ_σ₁ⱼ(i, j, k, grid, σ₁₁, σ₁₂) =
    (δxᶠᵃᵃ(i, j, k, grid, Δyᶜᶜ_σ, σ₁₁) + δyᵃᶜᵃ(i, j, k, grid, Δxᶠᶠ_σ, σ₁₂)) / Azᶠᶜᶜ(i, j, k, grid)

@inline old_∂ⱼ_σ₂ⱼ(i, j, k, grid, σ₂₂, σ₁₂) =
    (δxᶜᵃᵃ(i, j, k, grid, Δyᶠᶠ_σ, σ₁₂) + δyᵃᶠᵃ(i, j, k, grid, Δxᶜᶜ_σ, σ₂₂)) / Azᶜᶠᶜ(i, j, k, grid)

function set_smooth!(field, func; margin)
    grid = field.grid
    LX, LY, LZ = location(field)
    λ = λnodes(grid, LX(), LY(), LZ())
    φ = φnodes(grid, LX(), LY(), LZ())
    Nx, Ny, _ = size(grid)
    parent(field) .= 0
    for i in 1+margin:Nx-margin, j in 1+margin:Ny-margin
        field[i, j, 1] = func(λ[i], φ[j])
    end
    return nothing
end

function stress_power_budget(grid)
    rheology = ElastoViscoPlasticRheology()
    clock    = Clock(time=zero(eltype(grid)))

    u   = Field{Face,   Center, Nothing}(grid)
    v   = Field{Center, Face,   Nothing}(grid)
    σ₁₁ = Field{Center, Center, Nothing}(grid)
    σ₂₂ = Field{Center, Center, Nothing}(grid)
    σ₁₂ = Field{Face,   Face,   Nothing}(grid)

    fields = (; σ₁₁, σ₂₂, σ₁₂)

    Nx, Ny, Nz = size(grid)

    λ̂(λ) = (λ - 0)  / 60 * 2π   # longitude phase over the (0, 60) domain
    φ̂(φ) = (φ - 20) / 50 * 2π   # latitude  phase over the (20, 70) domain

    set_smooth!(u,   (λ, φ) -> sin(2λ̂(λ)) * cos(3φ̂(φ)); margin = 2)
    set_smooth!(v,   (λ, φ) -> cos(3λ̂(λ)) * sin(2φ̂(φ)); margin = 2)
    set_smooth!(σ₁₁, (λ, φ) -> sin(λ̂(λ))  * sin(2φ̂(φ)); margin = 2)
    set_smooth!(σ₂₂, (λ, φ) -> cos(2λ̂(λ)) * cos(φ̂(φ));  margin = 2)
    set_smooth!(σ₁₂, (λ, φ) -> sin(3λ̂(λ)) * cos(2φ̂(φ)); margin = 2)

    W_new = 0.0  # work by the metric-aware divergence
    W_old = 0.0  # work by the previous flux-form divergence
    D     = 0.0  # stress power σ:ε̇

    for i in 1:Nx, j in 1:Ny
        k = 1
        W_new += u[i, j, k] * ∂ⱼ_σ₁ⱼ(i, j, k, grid, rheology, clock, fields) * Azᶠᶜᶜ(i, j, k, grid)
        W_new += v[i, j, k] * ∂ⱼ_σ₂ⱼ(i, j, k, grid, rheology, clock, fields) * Azᶜᶠᶜ(i, j, k, grid)

        W_old += u[i, j, k] * old_∂ⱼ_σ₁ⱼ(i, j, k, grid, σ₁₁, σ₁₂) * Azᶠᶜᶜ(i, j, k, grid)
        W_old += v[i, j, k] * old_∂ⱼ_σ₂ⱼ(i, j, k, grid, σ₂₂, σ₁₂) * Azᶜᶠᶜ(i, j, k, grid)

        D += σ₁₁[i, j, k]     * strain_rate_xx(i, j, k, grid, u, v) * Azᶜᶜᶜ(i, j, k, grid)
        D += σ₂₂[i, j, k]     * strain_rate_yy(i, j, k, grid, u, v) * Azᶜᶜᶜ(i, j, k, grid)
        D += 2 * σ₁₂[i, j, k] * strain_rate_xy(i, j, k, grid, u, v) * Azᶠᶠᶜ(i, j, k, grid)
    end

    return W_new, W_old, D
end

relative_imbalance(W, D) = abs(W + D) / max(abs(W), abs(D))

@testset "Discrete energy budget of the stress divergence" begin
    # Evaluated at two resolutions on a curvilinear grid: the metric-aware operator is the
    # exact discrete adjoint of the strain-rate operator, so the energy identity is a pure
    # matrix-transpose identity that holds to machine precision regardless of resolution.
    # The previous flux-form operator is not the adjoint and carries a finite (percent-level)
    # imbalance at any usable resolution — orders of magnitude above roundoff.
    results = map((40, 80)) do N
        grid = LatitudeLongitudeGrid(size = (N, N),
                                     longitude = (0, 60),
                                     latitude = (20, 70),
                                     topology = (Bounded, Bounded, Flat),
                                     halo = (4, 4))

        W_new, W_old, D = stress_power_budget(grid)
        imbalance_new = relative_imbalance(W_new, D)
        imbalance_old = relative_imbalance(W_old, D)
        @info "energy budget" N imbalance_new imbalance_old
        (; new = imbalance_new, old = imbalance_old)
    end

    for r in results
        # Metric-aware formulation: exact discrete energy conservation (machine precision).
        @test r.new < 1e-10

        # Previous flux-form formulation: a finite energy imbalance survives.
        @test r.old > 1e-3

        # The metric-aware operator conserves energy by ~10⁶× more tightly at every resolution.
        @test r.new < 1e-6 * r.old
    end
end
