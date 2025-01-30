module Rheologies

export ViscousRheology, ElastoViscoPlasticRheology
export ∂ⱼ_σ₁ⱼ, ∂ⱼ_σ₂ⱼ, required_auxiliary_fields

using Oceananigans
using Oceananigans.Operators
using ClimaSeaIce: ice_mass

# Nothing rheology
required_auxiliary_fields(rheology, grid) = NamedTuple()
initialize_rheology!(model, rheology) = nothing
compute_stresses!(model, ice_dynamics, rheology, Δt) = nothing

@inline ∂ⱼ_σ₁ⱼ(i, j, k, grid, rheology, clock, fields) = zero(grid)
@inline ∂ⱼ_σ₂ⱼ(i, j, k, grid, rheology, clock, fields) = zero(grid)

# Nothing rheology or viscous rheology
@inline compute_time_stepᶠᶜᶜ(i, j, grid, Δt, rheology, substeps, fields) = Δt / substeps
@inline compute_time_stepᶜᶠᶜ(i, j, grid, Δt, rheology, substeps, fields) = Δt / substeps

# Fallback
@inline sum_of_forcing_u(i, j, k, grid, rheology, u_forcing, fields, Δt) = u_forcing(i, j, k, grid, fields)
@inline sum_of_forcing_v(i, j, k, grid, rheology, v_forcing, fields, Δt) = v_forcing(i, j, k, grid, fields)

include("viscous_rheology.jl")
include("strain_rates.jl")
include("elasto_visco_plastic_rheology.jl")

end