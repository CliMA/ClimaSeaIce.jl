module Rheologies

export ViscousRheology, ElastoViscoPlasticRheology
export ∂ⱼ_σ₁ⱼ, ∂ⱼ_σ₂ⱼ, Auxiliaries

using Adapt: Adapt
using Oceananigans: Oceananigans
using Oceananigans.Grids: AbstractGrid, Center, Face
using Oceananigans.Operators: Operators, Vᶜᶠᶜ, Vᶠᶜᶜ,
                              Axᶜᶜᶜ, Axᶠᶠᶜ, Ayᶜᶜᶜ, Ayᶠᶠᶜ, Azᶜᶜᶜ, Azᶠᶜᶜ, Azᶠᶠᶜ,
                              Δx_qᶜᶜᶜ, Δx_qᶜᶠᶜ, Δx_qᶠᶜᶜ, Δx_qᶠᶠᶜ,
                              Δy_qᶜᶜᶜ, Δy_qᶜᶠᶜ, Δy_qᶠᶜᶜ, Δy_qᶠᶠᶜ,
                              δxᶜᵃᵃ, δxᶜᶜᶜ, δxᶠᵃᵃ, δxᶠᶠᶜ,
                              δyᵃᶜᵃ, δyᵃᶠᵃ, δyᶜᶜᶜ, δyᶠᶠᶜ,
                              ℑxyᶜᶜᵃ, ℑxyᶠᶠᵃ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ
using Oceananigans.Utils: Utils, KernelParameters, configure_kernel

using ClimaSeaIce: ice_mass

struct Auxiliaries{F, K}
    fields :: F
    kernels :: K
end

# When adapted, only the fields need to be passed to the GPU.
# kernels operate only on the CPU.
Adapt.adapt_structure(to, a::Auxiliaries) =
    Auxiliaries(Adapt.adapt(to, a.fields), nothing)

"""
    Auxiliaries(rheology, grid)

A struct holding any auxiliary fields and kernels needed for the computation of
sea ice stresses.
"""
Auxiliaries(rheology, grid::AbstractGrid) = Auxiliaries(NamedTuple(), nothing)

import Oceananigans: prognostic_fields

# Nothing rheology
initialize_rheology!(model, rheology) = nothing
finalize_rheology!(fields, rheology) = nothing

compute_stresses!(dynamics, fields, grid, rheology, Δt) = nothing
prognostic_fields(mom, rheology) = NamedTuple()

# Nothing rheology or viscous rheology
@inline compute_substep_Δtᶠᶜᶜ(i, j, grid, Δt, rheology, substeps, fields) = Δt / substeps
@inline compute_substep_Δtᶜᶠᶜ(i, j, grid, Δt, rheology, substeps, fields) = Δt / substeps

# Fallback
@inline sum_of_forcing_u(i, j, k, grid, rheology, u_forcing, fields, Δt) = u_forcing(i, j, k, grid, fields)
@inline sum_of_forcing_v(i, j, k, grid, rheology, v_forcing, fields, Δt) = v_forcing(i, j, k, grid, fields)

@inline ice_stress_ux(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline ice_stress_uy(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline ice_stress_vx(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline ice_stress_vy(i, j, k, grid, ::Nothing, args...) = zero(grid)

include("ice_stress_divergence.jl")
include("viscous_rheology.jl")
include("elasto_visco_plastic_rheology.jl")

end