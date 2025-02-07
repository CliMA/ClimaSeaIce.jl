using Oceananigans.Grids: AbstractGrid, architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Architectures: convert_args
using Oceananigans.Utils: configure_kernel
using Oceananigans.TimeSteppers: store_field_tendencies!
using Oceananigans.ImmersedBoundaries: retrieve_surface_active_cells_map, mask_immersed_field_xy!

using ClimaSeaIce.Rheologies: _advance_bbm_stresses2!, _correct_bbm_stresses!
using Statistics: norm
mutable struct SplitExplicitSolver 
    substeps :: Int
end

"""
    SplitExplicitSolver(; substeps=120)

Creates a `SplitExplicitSolver` that controls the dynamical evolution of sea-ice momentum
by subcycling `substeps` times in between each thermodynamics / tracer advection time step.

The default number of substeps is 120.
"""
SplitExplicitSolver(; substeps=120) = SplitExplicitSolver(substeps)

const SplitExplicitMomentumEquation = SeaIceMomentumEquation{<:SplitExplicitSolver}

"""
    step_momentum!(model, rheology::AbstractExplicitRheology, Δt, χ)

function for stepping u and v in the case of _explicit_ solvers.
The sea-ice momentum equations are characterized by smaller time-scale than 
sea-ice thermodynamics and sea-ice tracer advection, therefore explicit rheologies require 
substepping over a set number of substeps.
"""
function step_momentum!(model, dynamics::SplitExplicitMomentumEquation, Δt, args...)

    grid = model.grid
    arch = architecture(grid)
    rheology = dynamics.rheology

    u, v = model.velocities
  
    ocean_velocities = dynamics.ocean_velocities
    clock = model.clock
    coriolis = dynamics.coriolis

    minimum_mass = dynamics.minimum_mass
    minimum_concentration = dynamics.minimum_concentration

    u_top_stress = dynamics.external_momentum_stresses.top.u
    v_top_stress = dynamics.external_momentum_stresses.top.v

    u_bottom_stress = dynamics.external_momentum_stresses.bottom.u
    v_bottom_stress = dynamics.external_momentum_stresses.bottom.v

    u_forcing = model.forcing.u
    v_forcing = model.forcing.v

    model_fields = merge(dynamics.auxiliary_fields, model.velocities, model.tracers,
                      (; h = model.ice_thickness, 
                         ℵ = model.ice_concentration, 
                         ρ = model.ice_density))

    u_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _u_velocity_step!)
    v_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _v_velocity_step!)

    substeps = dynamics.solver.substeps
    
    fill_halo_regions!(model.velocities)
    initialize_rheology!(model, dynamics.rheology)
    
    for substep in 1 : substeps
        # Compute stresses! depending on the particular rheology implementation

        not_converged = true
        iter = 1

        while not_converged 

            grid = model.grid
            arch = architecture(grid)

            ρᵢ   = model.ice_density
            d    = model.tracers.d
            u, v = model.velocities

            Nx, Ny, _ = size(grid)

            parameters = KernelParameters(-6:Nx+7, -6:Ny+7)

            # Pretty simple timestepping
            Δτ = Δt / substeps

            launch!(arch, grid, parameters, _advance_bbm_stresses2!, model_fields, grid, rheology, d, u, v, ρᵢ, Δτ)
            launch!(arch, grid, parameters, _correct_bbm_stresses!, model_fields, grid, rheology, d, ρᵢ, Δτ)
            # The momentum equations are solved using an alternating leap-frog algorithm
            # for u and v (used for the ocean - ice stresses and the coriolis term)
            # In even substeps we calculate uⁿ⁺¹ = f(vⁿ) and vⁿ⁺¹ = f(uⁿ⁺¹).
            # In odd substeps we switch and calculate vⁿ⁺¹ = f(uⁿ) and uⁿ⁺¹ = f(vⁿ⁺¹).
            if iseven(substep) 
                u_velocity_kernel!(u, grid, Δt, substeps, rheology, model_fields, 
                                   ocean_velocities, clock, coriolis,
                                   minimum_mass, minimum_concentration, 
                                   u_top_stress, u_bottom_stress, u_forcing)

                v_velocity_kernel!(v, grid, Δt, substeps, rheology, model_fields, 
                                   ocean_velocities, clock, coriolis, 
                                   minimum_mass, minimum_concentration,
                                   v_top_stress, v_bottom_stress, v_forcing)

            else
                v_velocity_kernel!(v, grid, Δt, substeps, rheology, model_fields, 
                                   ocean_velocities, clock, coriolis, 
                                   minimum_mass, minimum_concentration,
                                   v_top_stress, v_bottom_stress, v_forcing)


                u_velocity_kernel!(u, grid, Δt, substeps, rheology, model_fields, 
                                   ocean_velocities, clock, coriolis,
                                   minimum_mass, minimum_concentration, 
                                   u_top_stress, u_bottom_stress, u_forcing)
            end

            σ₁₁ₒ = interior(model_fields.σ₁₁ₒ)
            σ₂₂ₒ = interior(model_fields.σ₂₂ₒ)
            σ₁₂ₒ = interior(model_fields.σ₁₂ₒ)
            dₒ   = interior(model_fields.dₒ)

            σ₁₁ = interior(model_fields.σ₁₁)
            σ₂₂ = interior(model_fields.σ₂₂)
            σ₁₂ = interior(model_fields.σ₁₂)
            dᵢ  = interior(model_fields.d)

            nc₁ = norm(σ₁₁ₒ .- σ₁₁) > 1e-6
            nc₂ = norm(σ₂₂ₒ .- σ₂₂) > 1e-6
            nc₃ = norm(σ₁₂ₒ .- σ₁₂) > 1e-6            
            nc₄ = norm(dₒ   .- dᵢ ) > 1e-6

            # @show nc₁, nc₂, nc₃, nc₄

            not_converged = (nc₁ | nc₂ | nc₃ | nc₄) & (iter < 10000)

            iter += 1
            fill_my_halo_regions!(model.velocities)

            parent(model_fields.σ₁₁ₒ) .= parent(model_fields.σ₁₁)
            parent(model_fields.σ₂₂ₒ) .= parent(model_fields.σ₂₂)
            parent(model_fields.σ₁₂ₒ) .= parent(model_fields.σ₁₂)
            parent(model_fields.dₒ)   .= parent(model_fields.d)

            # TODO: This needs to be removed in some way!
        end

        @show iter

        parent(model_fields.uₙ)   .= parent(model_fields.u)
        parent(model_fields.vₙ)   .= parent(model_fields.v)
        parent(model_fields.dₙ)   .= parent(model_fields.d)
        parent(model_fields.σ₁₁ₙ) .= parent(model_fields.σ₁₁)
        parent(model_fields.σ₂₂ₙ) .= parent(model_fields.σ₂₂)
        parent(model_fields.σ₁₂ₙ) .= parent(model_fields.σ₁₂)
    end

    return nothing
end

@kernel function _fill_my_x_halo_regions!(u, v, Nx, Hx, Hy)
    j, _ = @index(Global, NTuple)
    j = j - Hy

    @inbounds u[1,    j, 1] = 0 # Impenetrability
    @inbounds u[Nx+1, j, 1] = 0 # Impenetrability
    @inbounds v[1,    j, 1] = 0 # Impenetrability
    @inbounds v[Nx+1, j, 1] = 0 # Impenetrability

    for i in 1:Hx-1
        @inbounds u[Nx+1+i, j, 1] = - u[Nx+1-i, j, 1]
        @inbounds v[Nx+1+i, j, 1] = - v[Nx+1-i, j, 1]
    end

    for i in 1:Hx
        @inbounds u[1-i, j, 1]  = - u[1+i, j, 1]
        @inbounds v[1-i, j, 1]  = - v[1+i, j, 1]
    end
end

@kernel function _fill_my_y_halo_regions!(u, v, Ny, Hx, Hy)
    i, _= @index(Global, NTuple)
    i = i - Hx
    @inbounds v[i, 1,    1] = 0 # Impenetrability
    @inbounds v[i, Ny+1, 1] = 0 # Impenetrability
    @inbounds u[i, 1,    1] = 0 # Impenetrability
    @inbounds u[i, Ny+1, 1] = 0 # Impenetrability

    for j in 1:Hy-1
        @inbounds v[i, Ny+1+j, 1] = - v[i, Ny+1-j, 1]
        @inbounds u[i, Ny+1+j, 1] = - u[i, Ny+1-j, 1]
    end
    
    for j in 1:Hy
        @inbounds v[i, 1-j, 1]  = - v[i, 1+j, 1]
        @inbounds u[i, 1-j, 1]  = - u[i, 1+j, 1]
    end
end

using Oceananigans.Grids: halo_size, architecture

@inline function fill_my_halo_regions!(velocities)
    u, v = velocities
    grid = u.grid
    arch = architecture(grid)
    Nx, Ny, _ = size(parent(u))
    nx, ny, _ = size(grid)
    Hx, Hy = halo_size(grid)

    launch!(arch, grid, (Ny, 1), _fill_my_x_halo_regions!, u, v, nx, Hx, Hy)
    launch!(arch, grid, (Nx, 1), _fill_my_y_halo_regions!, u, v, ny, Hx, Hy)
end

@kernel function _u_velocity_step!(u, grid, Δt, 
                                   substeps, rheology, 
                                   model_fields, ocean_velocities, 
                                   clock, coriolis, 
                                   minimum_mass, minimum_concentration,
                                   u_top_stress, u_bottom_stress, u_forcing)

    i, j = @index(Global, NTuple)

    mᵢ = ℑxyᶠᶠᵃ(i, j, 1, grid, ice_mass, model_fields.h, model_fields.ℵ, model_fields.ρ)
    ℵᵢ = ℑxyᶠᶠᵃ(i, j, 1, grid, model_fields.ℵ)

    Δτ = compute_time_stepᶠᶜᶜ(i, j, grid, Δt, rheology, substeps, model_fields) 
    Gu = u_velocity_tendency(i, j, grid, Δτ, rheology, model_fields, clock, coriolis, u_top_stress, u_bottom_stress, u_forcing)
   
    # Implicit part of the stress that depends linearly on the velocity
    τuᵢ = ( implicit_τx_coefficient(i, j, 1, grid, u_bottom_stress, clock, model_fields) 
          + implicit_τx_coefficient(i, j, 1, grid, u_top_stress, clock, model_fields)) / mᵢ * ℵᵢ 

    τuᵢ = ifelse(mᵢ ≤ 0, zero(grid), τuᵢ)
    uᴰ  = @inbounds (model_fields.uₙ[i, j, 1] + Δτ * Gu) / (1 + Δτ * τuᵢ) # dynamical velocity 
    uᶠ  = free_drift_u(i, j, 1, grid, ocean_velocities) # free drift velocity

    # If the ice mass or the ice concentration are below a certain threshold, 
    # the sea ice velocity is set to the free drift velocity
    sea_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)

    @inbounds u[i, j, 1] = ifelse(sea_ice, uᴰ, uᶠ)
end

@kernel function _v_velocity_step!(v, grid, Δt, 
                                   substeps, rheology, 
                                   model_fields, ocean_velocities, 
                                   clock, coriolis, 
                                   minimum_mass, minimum_concentration,
                                   v_top_stress, v_bottom_stress, v_forcing)

    i, j = @index(Global, NTuple)
    
    mᵢ = ℑxyᶠᶠᵃ(i, j, 1, grid, ice_mass, model_fields.h, model_fields.ℵ, model_fields.ρ)
    ℵᵢ = ℑxyᶠᶠᵃ(i, j, 1, grid, model_fields.ℵ)
    
    Δτ = compute_time_stepᶜᶠᶜ(i, j, grid, Δt, rheology, substeps, model_fields) 
    Gv = v_velocity_tendency(i, j, grid, Δτ, rheology, model_fields, clock, coriolis, v_top_stress, v_bottom_stress, v_forcing)

    # Implicit part of the stress that depends linearly on the velocity
    τvᵢ = ( implicit_τy_coefficient(i, j, 1, grid, v_bottom_stress, clock, model_fields)
          + implicit_τy_coefficient(i, j, 1, grid, v_top_stress, clock, model_fields)) / mᵢ * ℵᵢ 

    τvᵢ = ifelse(mᵢ ≤ 0, zero(grid), τvᵢ)

    vᴰ = @inbounds (model_fields.vₙ[i, j, 1] + Δτ * Gv) / (1 + Δτ * τvᵢ)# dynamical velocity 
    vᶠ = free_drift_v(i, j, 1, grid, ocean_velocities) # free drift velocity

    sea_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)

    @inbounds v[i, j, 1] = ifelse(sea_ice, vᴰ, vᶠ)
end