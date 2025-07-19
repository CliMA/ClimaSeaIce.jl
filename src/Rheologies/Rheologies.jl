module Rheologies

export ViscousRheology, ElastoViscoPlasticRheology, BrittleBinghamMaxwellRheology
export ∂ⱼ_σ₁ⱼ, ∂ⱼ_σ₂ⱼ, rheology_auxiliary_fields

using Oceananigans
using Oceananigans.Operators
using ClimaSeaIce: ice_mass

# Nothing rheology
initialize_rheology!(model, rheology) = nothing
compute_stresses!(model, dynamics, rheology, Δt, Ns) = nothing

# Nothing rheology or viscous rheology
@inline compute_substep_Δtᶠᶠᶜ(i, j, grid, Δt, rheology, substeps, fields) = Δt / substeps

# Fallback
@inline sum_of_forcing_u(i, j, k, grid, rheology, u_forcing, fields, Δt) = u_forcing(i, j, k, grid, fields)
@inline sum_of_forcing_v(i, j, k, grid, rheology, v_forcing, fields, Δt) = v_forcing(i, j, k, grid, fields)

# No needed extra tracers
rheology_prognostic_tracers(rheology) = ()

# No needed auxiliary fields
rheology_auxiliary_fields(rheology, grid) = NamedTuple()

@inline ice_stress_ux(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline ice_stress_uy(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline ice_stress_vx(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline ice_stress_vy(i, j, k, grid, ::Nothing, args...) = zero(grid)

include("ice_stress_divergence.jl")
include("viscous_rheology.jl")
include("elasto_visco_plastic_rheology.jl")
include("brittle_bingham_maxwell_rheology.jl")

end