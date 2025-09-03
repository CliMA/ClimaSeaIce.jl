module Rheologies

export ViscousRheology, ElastoViscoPlasticRheology
export ∂ⱼ_σ₁ⱼ, ∂ⱼ_σ₂ⱼ, required_auxiliaries

using Oceananigans
using Oceananigans.Operators
using Oceananigans.Grids: AbstractGrid
using ClimaSeaIce: ice_mass
using Adapt 

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

# Nothing rheology
required_auxiliaries(rheology, grid) = (; fields = NamedTuple())
initialize_rheology!(model, rheology) = nothing
compute_stresses!(dynamics, fields, grid, rheology, Δt) = nothing

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