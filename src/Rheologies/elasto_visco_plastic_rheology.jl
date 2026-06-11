using Adapt: Adapt
using KernelAbstractions: @kernel, @index
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.DistributedComputations: synchronize_communication!
using Oceananigans.Grids: AbstractGrid, halo_size

## The equations are solved in an iterative form following the EVP rheology of
## Kimmritz et al. (2017); doi: 10.1016/j.ocemod.2017.05.006
#
# Where:
# σᵢⱼ(u) = 2η ϵ̇ᵢⱼ + [(ζ - η) * (ϵ̇₁₁ + ϵ̇₂₂) - P / 2] δᵢⱼ
#
struct ElastoViscoPlasticRheology{FT, IP} <: AbstractRheology
    ice_compressive_strength :: FT # compressive strength
    ice_compaction_hardening :: FT # compaction hardening
    yield_curve_eccentricity :: FT # elliptic yield curve eccentricity
    minimum_plastic_stress :: FT # minimum plastic parameter (transitions to viscous behaviour)
    min_relaxation_parameter :: FT # minimum number of substeps expressed as the dynamic coefficient
    max_relaxation_parameter :: FT # maximum number of substeps expressed as the dynamic coefficient
    relaxation_strength :: FT # strength of the relaxation parameter
    pressure_formulation :: IP # formulation of ice pressure
    ElastoViscoPlasticRheology(P::FT, C::FT, e::FT, Δ_min::FT, α⁻::FT, α⁺::FT, c::FT, ip::IP) where {FT, IP}  =
        new{FT, IP}(P, C, e, Δ_min, α⁻, α⁺, c, ip)
end

function Base.show(io::IO, evpr::ElastoViscoPlasticRheology{FT}) where FT
    print(io, "ElastoViscoPlasticRheology{", FT, "}", '\n')
    print(io, "├── ice_compressive_strength: ", evpr.ice_compressive_strength, '\n')
    print(io, "├── ice_compaction_hardening: ", evpr.ice_compaction_hardening, '\n')
    print(io, "├── yield_curve_eccentricity: ", evpr.yield_curve_eccentricity, '\n')
    print(io, "├── minimum_plastic_stress: ", evpr.minimum_plastic_stress, '\n')
    print(io, "├── min_relaxation_parameter: ", evpr.min_relaxation_parameter, '\n')
    print(io, "├── max_relaxation_parameter: ", evpr.max_relaxation_parameter, '\n')
    print(io, "├── relaxation_strength: ", evpr.relaxation_strength, '\n')
    print(io, "└── pressure_formulation: ", summary(evpr.pressure_formulation))
end

struct ReplacementPressure end
struct IceStrength end

"""
    ElastoViscoPlasticRheology(FT = Oceananigans.defaults.FloatType;
                               ice_compressive_strength = 27500,
                               ice_compaction_hardening = 20,
                               yield_curve_eccentricity = 2,
                               minimum_plastic_stress = 2e-9,
                               min_relaxation_parameter = 50,
                               max_relaxation_parameter = 300,
                               relaxation_strength = π^2,
                               pressure_formulation = ReplacementPressure())

Construct an `ElastoViscoPlasticRheology` object representing a modified
elasto-visco-plastic rheology for sea-ice dynamics following
[Kimmritz et al. 2017](@cite Kimmritz2017).
The stress tensor is computed following the constitutive relation:
```math
σᵢⱼ = 2η ϵ̇ᵢⱼ + [(ζ - η) (ϵ̇₁₁ + ϵ̇₂₂) - P / 2] δᵢⱼ
```
where ``ϵ̇ᵢⱼ`` are the strain rates, ``η`` is the shear viscosity, ``ζ`` is the bulk viscosity,
and ``P`` is the ice strength (acting as the isotropic part of the stress tensor)
parameterized as ``P_\\star h \\exp[ - C ( 1 - ℵ )]`` where ``P_\\star`` is the `ice_compressive_strength`,
``C`` is the `ice_compaction_hardening`, ``h`` is the ice thickness, and ``ℵ`` is the ice concentration.

The stresses are substepped using a dynamic substepping coefficient ``α`` that is
spatially varying and computed dynamically as done by
[Kimmritz et al. 2017](@cite Kimmritz2017).
In particular: ``α = \\sqrt{γ²}``, where ``γ² = ζ c_α (Δt / mᵢ) / A_z`` is a stability parameter
with ``A_z`` is the area of the grid cell, ``mᵢ`` the ice mass, ``Δt`` the time step, and ``c_α`` a
numerical stability parameter which controls the strength of ``γ²``.

The stresses are substepped with:
```math
σᵢⱼᵖ⁺¹ = σᵢⱼᵖ + (σᵢⱼᵖ⁺¹ - σᵢⱼᵖ) / α
```

This formulation allows fast convergence in regions where ``α`` is small. Regions where
``α`` is large correspond to regions where the ice is more solid and the convergence is slower.
``α`` can be thought of as a "pseudo substep number" or a "relaxation parameter".
If we are using a subcycling solver, then if ``α`` ≪ number of substeps, the convergence is faster.

Arguments
=========

- `grid`: The computational grid.

Keyword Arguments
=================

- `ice_compressive_strength`: Parameter expressing compressive strength
                              (in N m⁻²). Default: `27500`.
- `ice_compaction_hardening`: Exponent coefficient for compaction hardening.
                              Default: `20`.
- `yield_curve_eccentricity`: Eccentricity of the elliptic yield curve.
                              Default: `2`.
- `minimum_plastic_stress`: Minimum value for the visco-plastic parameter.
                            Limits the maximum viscosity of the ice, transitioning
                            the rheology from plastic to viscous behavior.
                            Default: `2e-9`.
- `min_relaxation_parameter`: Minimum value for the relaxation parameter `α`.
                              Default: `50`.
- `max_relaxation_parameter`: Maximum value for the relaxation parameter `α`.
                              Default: `300`.
- `relaxation_strength`: Parameter controlling the strength of the relaxation
                         parameter. The maximum value is `π²`; see
                         [Kimmritz et al. 2017](@cite Kimmritz2017).
                         Default: `π²`.
- `pressure_formulation`: Either `ReplacementPressure()` or `IceStrength()`.
                          The replacement-pressure formulation avoids ice motion
                          in the absence of forcing. Default:
                          `ReplacementPressure()`.

References
==========

- Kimmritz, M., Losch, M., and Danilov, S. (2017). A comparison of viscous-plastic sea ice solvers with and without
  replacement pressure. Ocean Modelling, 115, 59-69. doi:10.1016/j.ocemod.2017.05.006.
"""
function ElastoViscoPlasticRheology(FT::DataType = Oceananigans.defaults.FloatType;
                                    ice_compressive_strength = 27500,
                                    ice_compaction_hardening = 20,
                                    yield_curve_eccentricity = 2,
                                    minimum_plastic_stress = 2e-9,
                                    min_relaxation_parameter = 50,
                                    max_relaxation_parameter = 300,
                                    relaxation_strength = π^2,
                                    pressure_formulation = ReplacementPressure())

    return ElastoViscoPlasticRheology(convert(FT, ice_compressive_strength),
                                      convert(FT, ice_compaction_hardening),
                                      convert(FT, yield_curve_eccentricity),
                                      convert(FT, minimum_plastic_stress),
                                      convert(FT, min_relaxation_parameter),
                                      convert(FT, max_relaxation_parameter),
                                      convert(FT, relaxation_strength),
                                      pressure_formulation)
end

# Extend Auxiliaries to hold auxiliaries for the ElastoViscoPlasticRheology
function Auxiliaries(r::ElastoViscoPlasticRheology, grid::AbstractGrid)

    arch       = architecture(grid)
    Nx, Ny, _  = size(grid)
    Hx, Hy, _  = halo_size(grid)
    parameters = KernelParameters(-Hx+2:Nx+Hx-1, -Hy+2:Ny+Hy-1)

    σ₁₁ = Field{Center, Center, Nothing}(grid)
    σ₂₂ = Field{Center, Center, Nothing}(grid)
    σ₁₂ = Field{Face,   Face,   Nothing}(grid)
    uⁿ  = Field{Face,   Center, Nothing}(grid)
    vⁿ  = Field{Center, Face,   Nothing}(grid)
    P   = Field{Center, Center, Nothing}(grid)
    α   = Field{Center, Center, Nothing}(grid) # Dynamic substeps a la Kimmritz et al. (2017)
    Δ   = Field{Center, Center, Nothing}(grid)

    # Viscosities
    ζᶠᶠᶜ = Field{Face,   Face,   Nothing}(grid)
    ζᶜᶜᶜ = Field{Center, Center, Nothing}(grid)

    # An initial (safe) educated guess
    fill!(α, r.max_relaxation_parameter)

    _viscosity_kernel! = configure_kernel(arch, grid, parameters, _compute_evp_viscosities!)[1]
    _stresses_kernel!  = configure_kernel(arch, grid, parameters, _compute_evp_stresses!)[1]

    parameters = KernelParameters(size(P.data)[1:2], P.data.offsets[1:2])
    _initialize_rhology! = configure_kernel(arch, grid, parameters, _initialize_evp_rhology!)[1]

    fields  = (; σ₁₁, σ₂₂, σ₁₂, ζᶠᶠᶜ, ζᶜᶜᶜ, Δ, α, uⁿ, vⁿ, P)
    kernels = (; _viscosity_kernel!, _stresses_kernel!, _initialize_rhology!)

    return Auxiliaries(fields, kernels)
end

# Extend the `adapt_structure` function for the ElastoViscoPlasticRheology
Adapt.adapt_structure(to, r::ElastoViscoPlasticRheology) =
    ElastoViscoPlasticRheology(Adapt.adapt(to, r.ice_compressive_strength),
                               Adapt.adapt(to, r.ice_compaction_hardening),
                               Adapt.adapt(to, r.yield_curve_eccentricity),
                               Adapt.adapt(to, r.minimum_plastic_stress),
                               Adapt.adapt(to, r.min_relaxation_parameter),
                               Adapt.adapt(to, r.max_relaxation_parameter),
                               Adapt.adapt(to, r.relaxation_strength),
                               Adapt.adapt(to, r.pressure_formulation))

"""
    initialize_rheology!(model, rheology::ElastoViscoPlasticRheology)

Initialize the elasto-visco-plastic rheology.
In this step we calculate the ice strength given the ice mass (thickness and concentration).
"""
function initialize_rheology!(model, rheology::ElastoViscoPlasticRheology)
    h = model.ice_thickness
    ℵ = model.ice_concentration

    P★ = rheology.ice_compressive_strength
    C  = rheology.ice_compaction_hardening

    u, v    = model.velocities
    fields  = model.dynamics.auxiliaries.fields
    kernels = model.dynamics.auxiliaries.kernels
    kernels._initialize_rhology!(fields, model.grid, P★, C, h, ℵ, u, v)

    synchronize_communication!(fields.σ₁₁)
    synchronize_communication!(fields.σ₁₂)
    synchronize_communication!(fields.σ₂₂)

    return nothing
end

@kernel function _initialize_evp_rhology!(fields, grid, P★, C, h, ℵ, u, v)
    i, j = @index(Global, NTuple)
    @inbounds fields.P[i, j, 1]  = ice_strength(i, j, 1, grid, P★, C, h, ℵ)
    @inbounds fields.uⁿ[i, j, 1] = u[i, j, 1]
    @inbounds fields.vⁿ[i, j, 1] = v[i, j, 1]
end

# The parameterization for an `ElastoViscoPlasticRheology`
@inline ice_strength(i, j, k, grid, P★, C, h, ℵ) = @inbounds P★ * h[i, j, k] * exp(- C * (1 - ℵ[i, j, k]))

# Specific compute stresses for the EVP rheology
function compute_stresses!(dynamics, fields, grid, rheology::ElastoViscoPlasticRheology, Δt)

    h  = fields.h
    ρᵢ = fields.ρ
    ℵ  = fields.ℵ
    u  = fields.u
    v  = fields.v

    dynamics.auxiliaries.kernels._viscosity_kernel!(fields, grid, rheology, u, v)
    dynamics.auxiliaries.kernels._stresses_kernel!(fields, grid, rheology, u, v, h, ℵ, ρᵢ, Δt)

    return nothing
end

@inline strain_rate_xx(i, j, k, grid, u, v) =  δxᶜᵃᵃ(i, j, k, grid, Δy_qᶠᶜᶜ, u) / Azᶜᶜᶜ(i, j, k, grid)
@inline strain_rate_yy(i, j, k, grid, u, v) =  δyᵃᶜᵃ(i, j, k, grid, Δx_qᶜᶠᶜ, v) / Azᶜᶜᶜ(i, j, k, grid)
@inline strain_rate_xy(i, j, k, grid, u, v) = (δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶠᶜ, v) + δyᵃᶠᵃ(i, j, k, grid, Δx_qᶠᶜᶜ, u)) / Azᶠᶠᶜ(i, j, k, grid) / 2

@kernel function _compute_evp_viscosities!(fields, grid, rheology, u, v)
    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    Δm  = rheology.minimum_plastic_stress

    # Extract auxiliary fields
    P = fields.P

    # Strain rates
    ϵ̇₁₁ᶜᶜᶜ = strain_rate_xx(i, j, kᴺ, grid, u, v)
    ϵ̇₂₂ᶜᶜᶜ = strain_rate_yy(i, j, kᴺ, grid, u, v)
    ϵ̇₁₂ᶠᶠᶜ = strain_rate_xy(i, j, kᴺ, grid, u, v)
    ϵ̇₁₁ᶠᶠᶜ = ℑxyᶠᶠᵃ(i, j, kᴺ, grid, strain_rate_xx, u, v)
    ϵ̇₂₂ᶠᶠᶜ = ℑxyᶠᶠᵃ(i, j, kᴺ, grid, strain_rate_yy, u, v)
    ϵ̇₁₂ᶜᶜᶜ = ℑxyᶜᶜᵃ(i, j, kᴺ, grid, strain_rate_xy, u, v)

    # Ice divergence
    δᶜᶜᶜ = ϵ̇₁₁ᶜᶜᶜ + ϵ̇₂₂ᶜᶜᶜ
    δᶠᶠᶜ = ϵ̇₁₁ᶠᶠᶜ + ϵ̇₂₂ᶠᶠᶜ

    # Ice shear (at Centers)
    sᶜᶜᶜ = sqrt((ϵ̇₁₁ᶜᶜᶜ - ϵ̇₂₂ᶜᶜᶜ)^2 + 4ϵ̇₁₂ᶜᶜᶜ^2)
    sᶠᶠᶜ = sqrt((ϵ̇₁₁ᶠᶠᶜ - ϵ̇₂₂ᶠᶠᶜ)^2 + 4ϵ̇₁₂ᶠᶠᶜ^2)

    # Visco - Plastic parameter
    # if Δ is very small we assume a linear viscous response
    # adding a minimum Δ_min (at Centers)
    Δᶜᶜᶜ = max(sqrt(δᶜᶜᶜ^2 + sᶜᶜᶜ^2 * e⁻²), Δm)
    Δᶠᶠᶜ = max(sqrt(δᶠᶠᶜ^2 + sᶠᶠᶜ^2 * e⁻²), Δm)
    Pᶜᶜᶜ = @inbounds P[i, j, 1]
    Pᶠᶠᶜ = ℑxyᶠᶠᵃ(i, j, 1, grid, P)

    @inbounds fields.ζᶠᶠᶜ[i, j, 1] = Pᶠᶠᶜ / 2Δᶠᶠᶜ
    @inbounds fields.ζᶜᶜᶜ[i, j, 1] = Pᶜᶜᶜ / 2Δᶜᶜᶜ
    @inbounds fields.Δ[i, j, 1] = Δᶜᶜᶜ
end

function finalize_rheology!(fields, ::ElastoViscoPlasticRheology)
    fill_halo_regions!(fields.σ₁₁; async=true)
    fill_halo_regions!(fields.σ₁₂; async=true)
    fill_halo_regions!(fields.σ₂₂; async=true)
    return nothing
end

@inline ice_pressure(i, j, k, grid, ::IceStrength, r, fields) = @inbounds fields.P[i, j, k]

@inline function ice_pressure(i, j, k, grid, ::ReplacementPressure, r, fields)
    Pᶜᶜᶜ = @inbounds fields.P[i, j, k]
    Δᶜᶜᶜ = @inbounds fields.Δ[i, j, k]
    Δm   = r.minimum_plastic_stress
    return Pᶜᶜᶜ * Δᶜᶜᶜ / (Δᶜᶜᶜ + Δm)
end

# Compute the visco-plastic stresses for a slab sea ice model.
# The function updates the internal stress variables `σ₁₁`, `σ₂₂`, and `σ₁₂` in the `rheology` object
# following the αEVP formulation of Kimmritz et al. (2017).
@kernel function _compute_evp_stresses!(fields, grid, rheology, u, v, h, ℵ, ρᵢ, Δt)
    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    e⁻² = rheology.yield_curve_eccentricity^(-2)
    α⁺  = rheology.max_relaxation_parameter
    α⁻  = rheology.min_relaxation_parameter
    cα  = rheology.relaxation_strength
    ip  = rheology.pressure_formulation

    σ₁₁ = fields.σ₁₁
    σ₂₂ = fields.σ₂₂
    σ₁₂ = fields.σ₁₂
    α   = fields.α

    # Strain rates
    ϵ̇₁₁ = strain_rate_xx(i, j, kᴺ, grid, u, v)
    ϵ̇₂₂ = strain_rate_yy(i, j, kᴺ, grid, u, v)
    ϵ̇₁₂ = strain_rate_xy(i, j, kᴺ, grid, u, v)

    ζᶜᶜᶜ = @inbounds fields.ζᶜᶜᶜ[i, j, 1]
    ζᶠᶠᶜ = @inbounds fields.ζᶠᶠᶜ[i, j, 1]

    # replacement pressure?
    Pᵣ = ice_pressure(i, j, 1, grid, ip, rheology, fields)

    ηᶜᶜᶜ = ζᶜᶜᶜ * e⁻²
    ηᶠᶠᶜ = ζᶠᶠᶜ * e⁻²

    # σ(uᵖ): the tangential stress depends only shear viscosity
    # while the compressive stresses depend on the bulk viscosity and the ice strength
    σ₁₁ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₁₁ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2)
    σ₂₂ᵖ⁺¹ = 2 * ηᶜᶜᶜ * ϵ̇₂₂ + ((ζᶜᶜᶜ - ηᶜᶜᶜ) * (ϵ̇₁₁ + ϵ̇₂₂) - Pᵣ / 2)
    σ₁₂ᵖ⁺¹ = 2 * ηᶠᶠᶜ * ϵ̇₁₂

    mᵢᶜᶜᶜ = ice_mass(i, j, 1, grid, h, ℵ, ρᵢ)
    mᵢᶠᶠᶜ = ℑxyᶠᶠᵃ(i, j, 1, grid, ice_mass, h, ℵ, ρᵢ)

    # Update coefficients for substepping using dynamic substepping
    # with spatially varying coefficients as done by Kimmritz et al. (2017)
    γ²ᶜᶜᶜ = ζᶜᶜᶜ * cα * Δt / mᵢᶜᶜᶜ / Azᶜᶜᶜ(i, j, 1, grid)
    γ²ᶜᶜᶜ = ifelse(isnan(γ²ᶜᶜᶜ), α⁺^2, γ²ᶜᶜᶜ) # In case both ζᶜᶜᶜ and mᵢᶜᶜᶜ are zero
    γᶜᶜᶜ  = clamp(sqrt(γ²ᶜᶜᶜ), α⁻, α⁺)

    γ²ᶠᶠᶜ = ζᶠᶠᶜ * cα * Δt / mᵢᶠᶠᶜ / Azᶠᶠᶜ(i, j, 1, grid)
    γ²ᶠᶠᶜ = ifelse(isnan(γ²ᶠᶠᶜ), α⁺^2, γ²ᶠᶠᶜ) # In case both ζᶠᶠᶜ and mᵢᶠᶠᶜ are zero
    γᶠᶠᶜ  = clamp(sqrt(γ²ᶠᶠᶜ), α⁻, α⁺)

    @inbounds begin
        # Compute the new stresses and store the value of the
        # dynamic substepping coefficient α
        σ₁₁★ = (σ₁₁ᵖ⁺¹ - σ₁₁[i, j, 1]) / γᶜᶜᶜ
        σ₂₂★ = (σ₂₂ᵖ⁺¹ - σ₂₂[i, j, 1]) / γᶜᶜᶜ
        σ₁₂★ = (σ₁₂ᵖ⁺¹ - σ₁₂[i, j, 1]) / γᶠᶠᶜ

        σ₁₁[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, σ₁₁★, zero(grid))
        σ₂₂[i, j, 1] += ifelse(mᵢᶜᶜᶜ > 0, σ₂₂★, zero(grid))
        σ₁₂[i, j, 1] += ifelse(mᵢᶠᶠᶜ > 0, σ₁₂★, zero(grid))
          α[i, j, 1]  = γᶜᶜᶜ
    end
end

#####
##### Internal stress divergence for the EVP model
#####

# Here we extend all the functions that a rheology model needs to support:
@inline ice_stress_ux(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = @inbounds fields.σ₁₁[i, j, k]
@inline ice_stress_vx(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = @inbounds fields.σ₁₂[i, j, k]
@inline ice_stress_uy(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = @inbounds fields.σ₁₂[i, j, k]
@inline ice_stress_vy(i, j, k, grid, ::ElastoViscoPlasticRheology, clock, fields) = @inbounds fields.σ₂₂[i, j, k]

# To help convergence to the right velocities
@inline compute_substep_Δtᶠᶜᶜ(i, j, grid, Δt, ::ElastoViscoPlasticRheology, substeps, fields) = Δt / ℑxᶠᵃᵃ(i, j, 1, grid, fields.α)
@inline compute_substep_Δtᶜᶠᶜ(i, j, grid, Δt, ::ElastoViscoPlasticRheology, substeps, fields) = Δt / ℑyᵃᶠᵃ(i, j, 1, grid, fields.α)

#####
##### Numerical forcing to help convergence
#####

@inline function sum_of_forcing_u(i, j, k, grid, ::ElastoViscoPlasticRheology, u_forcing, fields, Δt)
    user_forcing = u_forcing(i, j, k, grid, fields)
    rheology_forcing = @inbounds (fields.uⁿ[i, j, k] - fields.u[i, j, k]) / Δt / ℑxᶠᵃᵃ(i, j, k, grid, fields.α)
    return user_forcing + rheology_forcing
end

@inline function sum_of_forcing_v(i, j, k, grid, ::ElastoViscoPlasticRheology, v_forcing, fields, Δt)
    user_forcing = v_forcing(i, j, k, grid, fields)
    rheology_forcing = @inbounds (fields.vⁿ[i, j, k] - fields.v[i, j, k]) / Δt / ℑyᵃᶠᵃ(i, j, k, grid, fields.α)
    return user_forcing + rheology_forcing
end
