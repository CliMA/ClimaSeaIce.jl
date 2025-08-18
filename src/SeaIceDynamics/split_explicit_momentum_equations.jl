using Oceananigans.Grids: AbstractGrid, architecture, halo_size
using Oceananigans.BoundaryConditions: fill_halo_regions!, fill_halo_size, fill_halo_offset
using Oceananigans.Utils: configure_kernel
using Oceananigans.Fields: instantiated_location, boundary_conditions
using Oceananigans.ImmersedBoundaries: peripheral_node

struct SplitExplicitSolver 
    substeps :: Int
end

"""
    SplitExplicitSolver(; substeps=120)

Creates a `SplitExplicitSolver` that controls the dynamical evolution of sea-ice momentum
by subcycling `substeps` times in between each ice_thermodynamics / tracer advection time step.

The default number of substeps is 120.
"""
SplitExplicitSolver(; substeps=120) = SplitExplicitSolver(substeps)

const SplitExplicitMomentumEquation = SeaIceMomentumEquation{<:SplitExplicitSolver}

"""
    time_step_momentum!(model, rheology::AbstractExplicitRheology, Δt)

function for stepping u and v in the case of _explicit_ solvers.
The sea-ice momentum equations are characterized by smaller time-scale than 
sea-ice ice_thermodynamics and sea-ice tracer advection, therefore explicit rheologies require 
substepping over a set number of substeps.
"""
function time_step_momentum!(model, dynamics::SplitExplicitMomentumEquation, Δt)

    grid = model.grid
    arch = architecture(grid)
    rheology = dynamics.rheology
    Nx, Ny, Nz = size(grid)
    Hx, Hy, _  = halo_size(grid)

    # Params for filling halo regions
    params_x = KernelParameters(-Hy:Ny+Hy, Nz:Nz)
    params_y = KernelParameters(-Hx:Nx+Hx, Nz:Nz)
    
    u, v = model.velocities
  
    free_drift = dynamics.free_drift
    clock = model.clock
    coriolis = dynamics.coriolis

    minimum_mass = dynamics.minimum_mass
    minimum_concentration = dynamics.minimum_concentration

    top_stress = dynamics.external_momentum_stresses.top
    bottom_stress = dynamics.external_momentum_stresses.bottom

    u_forcing = model.forcing.u
    v_forcing = model.forcing.v
    u_immersed_bc = u.boundary_conditions.immersed
    v_immersed_bc = v.boundary_conditions.immersed

    model_fields = merge(dynamics.auxiliary_fields, model.velocities, 
                      (; h = model.ice_thickness, 
                         ℵ = model.ice_concentration, 
                         ρ = model.ice_density))

    active_cells_map = Oceananigans.Grids.get_active_column_map(grid)

    u_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _u_velocity_step!; active_cells_map)
    v_velocity_kernel!, _ = configure_kernel(arch, grid, :xy, _v_velocity_step!; active_cells_map)

    substeps = dynamics.solver.substeps
    initialize_rheology!(model, dynamics.rheology)

    u_args = (u, grid, Δt, substeps, rheology, model_fields, 
              free_drift, clock, coriolis,
              minimum_mass, minimum_concentration, 
              u_immersed_bc, top_stress, bottom_stress, u_forcing)

    v_args = (v, grid, Δt, substeps, rheology, model_fields, 
              free_drift, clock, coriolis, 
              minimum_mass, minimum_concentration,
              v_immersed_bc, top_stress, bottom_stress, v_forcing)

    u_fill_halo_args = (u.data, u.boundary_conditions, u.indices, (Face(), Center(), nothing), grid)
    v_fill_halo_args = (v.data, v.boundary_conditions, v.indices, (Center(), Face(), nothing), grid)
    stresses_args    = (model_fields, dynamics, grid, rheology, Δt)

    GC.@preserve v_args u_args u_fill_halo_args v_fill_halo_args stresses_args begin
        # We need to perform ~150 time-steps which means
        # launching ~300 very small kernels: we are limited by
        # latency of argument conversion to GPU-compatible values.
        # To alleviate this penalty we convert first and then we substep!
        converted_u_args = Oceananigans.Architectures.convert_to_device(arch, u_args)
        converted_v_args = Oceananigans.Architectures.convert_to_device(arch, v_args)
        converted_u_halo = Oceananigans.Architectures.convert_to_device(arch, u_fill_halo_args)
        converted_v_halo = Oceananigans.Architectures.convert_to_device(arch, v_fill_halo_args)
        converted_stresses_args = Oceananigans.Architectures.convert_to_device(arch, stresses_args)

        for substep in 1 : substeps
            # Compute stresses! depending on the particular rheology implementation
            compute_stresses!(converted_stresses_args...)

            # The momentum equations are solved using an alternating leap-frog algorithm
            # for u and v (used for the ocean - ice stresses and the coriolis term)
            # In even substeps we calculate uⁿ⁺¹ = f(vⁿ) and vⁿ⁺¹ = f(uⁿ⁺¹).
            # In odd substeps we switch and calculate vⁿ⁺¹ = f(uⁿ) and uⁿ⁺¹ = f(vⁿ⁺¹).
            if iseven(substep) 
                u_velocity_kernel!(converted_u_args...)
                v_velocity_kernel!(converted_v_args...)
            else
                v_velocity_kernel!(converted_v_args...)
                u_velocity_kernel!(converted_u_args...)
            end

            fill_halo_regions!(converted_u_halo...)
            fill_halo_regions!(converted_v_halo...)
        end
    end

    return nothing
end

@kernel function _u_velocity_step!(u, grid, Δt, 
                                   substeps, rheology, 
                                   model_fields, free_drift, 
                                   clock, coriolis, 
                                   minimum_mass, minimum_concentration,
                                   u_immersed_bc, u_top_stress, u_bottom_stress, u_forcing)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) 

    mᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, ice_mass, model_fields.h, model_fields.ℵ, model_fields.ρ)
    ℵᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, model_fields.ℵ)

    Δτ = compute_substep_Δtᶠᶜᶜ(i, j, grid, Δt, rheology, substeps, model_fields) 
    Gu = u_velocity_tendency(i, j, grid, Δτ, rheology, model_fields, clock, coriolis, u_immersed_bc, u_top_stress, u_bottom_stress, u_forcing)
   
    # Implicit part of the stress that depends linearly on the velocity
    τuᵢ = ( implicit_τx_coefficient(i, j, kᴺ, grid, u_bottom_stress, clock, model_fields) 
          - implicit_τx_coefficient(i, j, kᴺ, grid, u_top_stress, clock, model_fields)) / mᵢ * ℵᵢ 

    τuᵢ = ifelse(mᵢ ≤ 0, zero(grid), τuᵢ)
    uᴰ  = @inbounds (u[i, j, 1] + Δτ * Gu) / (1 + Δτ * τuᵢ) # dynamical velocity 
    uᶠ  = free_drift_u(i, j, kᴺ, grid, free_drift, clock, model_fields) # free drift velocity

    # If the ice mass or the ice concentration are below a certain threshold, 
    # the sea ice velocity is set to the free drift velocity
    sea_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)
    active  = !peripheral_node(i, j, kᴺ, grid, Face(), Center(), Center())

    @inbounds u[i, j, 1] = ifelse(sea_ice, uᴰ, uᶠ) * active
end

@kernel function _v_velocity_step!(v, grid, Δt, 
                                   substeps, rheology, 
                                   model_fields, free_drift, 
                                   clock, coriolis, 
                                   minimum_mass, minimum_concentration,
                                   v_immersed_bc, v_top_stress, v_bottom_stress, v_forcing)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3) 

    mᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, ice_mass, model_fields.h, model_fields.ℵ, model_fields.ρ)
    ℵᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, model_fields.ℵ)
    
    Δτ = compute_substep_Δtᶜᶠᶜ(i, j, grid, Δt, rheology, substeps, model_fields) 
    Gv = v_velocity_tendency(i, j, grid, Δτ, rheology, model_fields, clock, coriolis, v_immersed_bc, v_top_stress, v_bottom_stress, v_forcing)

    # Implicit part of the stress that depends linearly on the velocity
    τvᵢ = ( implicit_τy_coefficient(i, j, kᴺ, grid, v_bottom_stress, clock, model_fields)
          - implicit_τy_coefficient(i, j, kᴺ, grid, v_top_stress, clock, model_fields)) / mᵢ * ℵᵢ 

    τvᵢ = ifelse(mᵢ ≤ 0, zero(grid), τvᵢ)

    vᴰ = @inbounds (v[i, j, 1] + Δτ * Gv) / (1 + Δτ * τvᵢ)# dynamical velocity 
    vᶠ = free_drift_v(i, j, kᴺ, grid, free_drift, clock, model_fields)  # free drift velocity

    # If the ice mass or the ice concentration are below a certain threshold, 
    # the sea ice velocity is set to the free drift velocity
    sea_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)
    active  = !peripheral_node(i, j, kᴺ, grid, Center(), Face(), Center())

    @inbounds v[i, j, 1] = ifelse(sea_ice, vᴰ, vᶠ) * active
end
