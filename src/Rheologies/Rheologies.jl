module Rheologies

export ViscousRheology, ElastoViscoPlasticRheology
export FreeSlip, NoSlip
export ∂ⱼ_σ₁ⱼ, ∂ⱼ_σ₂ⱼ, Auxiliaries

using Adapt: Adapt
using Oceananigans: Oceananigans
using Oceananigans.Fields: Field
using Oceananigans.Grids: AbstractGrid, Center, Face
using Oceananigans.Operators: Azᶜᶜᶜ, Azᶜᶠᶜ, Azᶠᶜᶜ, Azᶠᶠᶜ,
                              Δx_qᶜᶠᶜ, Δy_qᶠᶜᶜ,
                              Δxᶜᶜᶜ, Δxᶜᶠᶜ, Δxᶠᶜᶜ, Δxᶠᶠᶜ,
                              Δyᶜᶜᶜ, Δyᶜᶠᶜ, Δyᶠᶜᶜ, Δyᶠᶠᶜ,
                              δxᶜᵃᵃ, δxᶜᶜᶜ, δxᶠᵃᵃ, δxᶠᶠᶜ,
                              δyᵃᶜᵃ, δyᵃᶠᵃ, δyᶜᶜᶜ, δyᶠᶠᶜ,
                              ℑxyᶜᶜᵃ, ℑxyᶠᶠᵃ, ℑxᶠᵃᵃ, ℑyᵃᶠᵃ
using Oceananigans.Utils: KernelParameters, configure_kernel

using ..ClimaSeaIce: ice_mass

abstract type AbstractRheology end

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
initialize_rheology!(model, rheology) = nothing
finalize_rheology!(fields, rheology) = nothing

compute_stresses!(kernels, fields, grid, rheology, Δt, u_immersed_bc, v_immersed_bc) = nothing

# Rheologies that do not store stresses in auxiliary fields have no stress kernels to restrict
stress_kernel_ranges(rheology, grid) = nothing
mapped_stress_kernels(kernels, rheology, arch, grid, active_cells_map) = kernels

# A kernel launched over an empty map does nothing
struct NoKernel end
@inline (::NoKernel)(args...) = nothing

"""
    configure_mapped_kernel(arch, grid, kernel!, active_cells_map)

Configure `kernel!` to run only over the `(i, j)` indices listed in `active_cells_map`.
Returns a `NoKernel` (which does nothing when called) if `active_cells_map` is empty.
"""
function configure_mapped_kernel(arch, grid, kernel!, active_cells_map)
    isempty(active_cells_map) && return NoKernel()
    return first(configure_kernel(arch, grid, :xy, kernel!; active_cells_map))
end
Oceananigans.prognostic_fields(mom, ::AbstractRheology) = NamedTuple()

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

include("lateral_boundary_conditions.jl")
include("ice_stress_divergence.jl")
include("viscous_rheology.jl")
include("elasto_visco_plastic_rheology.jl")

end
