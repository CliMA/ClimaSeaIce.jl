using Oceananigans: instantiated_location
using Oceananigans.Architectures: convert_to_device, architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.DistributedComputations: DistributedGrid
using Oceananigans.Grids: AbstractGrid, halo_size, halo_size, topology, with_halo, peripheral_node,
                          LeftConnected, RightConnected, FullyConnected,
                          RightCenterFolded, RightFaceFolded,
                          LeftConnectedRightCenterFolded, LeftConnectedRightFaceFolded,
                          LeftConnectedRightCenterConnected, LeftConnectedRightFaceConnected
using Oceananigans.Models.HydrostaticFreeSurfaceModels.SplitExplicitFreeSurfaces: split_explicit_kernel_size
using Oceananigans.Utils: configure_kernel

const ConnectedTopology = Union{LeftConnected, RightConnected, FullyConnected,
                                RightCenterFolded, RightFaceFolded,
                                LeftConnectedRightCenterFolded, LeftConnectedRightFaceFolded,
                                LeftConnectedRightCenterConnected, LeftConnectedRightFaceConnected}

struct SplitExplicitSolver{I, K}
    substeps :: I
    kernel_parameters :: K
end

"""
    SplitExplicitSolver(grid::AbstractGrid; substeps=120)

Creates a `SplitExplicitSolver` that controls the dynamical evolution of sea-ice momentum
by subcycling `substeps` times in between each ice_thermodynamics / tracer advection time step.

The default number of substeps is 120.
"""
SplitExplicitSolver(grid::AbstractGrid; substeps=120) = SplitExplicitSolver(substeps, :xy)

# When no grid is provided, we assume a serial grid with default kernel parameters
SplitExplicitSolver(; substeps=120) = SplitExplicitSolver(substeps, :xy)

const SplitExplicitMomentumEquation = SeaIceMomentumEquation{<:SplitExplicitSolver}

# Shenanigans for extending the halos in Distributed grids

function SplitExplicitSolver(grid::DistributedGrid; substeps=120)
    Nx, Ny, _ = size(grid)
    Hx, Hy, _ = halo_size(grid)
    TX, TY, _ = topology(grid)
    kernel_sizes = map(split_explicit_kernel_size, (TX, TY), (Nx, Ny), (Hx, Hy))
    return SplitExplicitSolver(substeps, KernelParameters(kernel_sizes...))
end

function maybe_extended_grid(mom::SplitExplicitMomentumEquation, grid::DistributedGrid)
    old_halos = halo_size(grid)
    Nsubsteps = mom.solver.substeps
    TX, TY, _ = topology(grid)
    Hx = TX() isa ConnectedTopology ? max(2Nsubsteps + 3, old_halos[1]) : old_halos[1]
    Hy = TY() isa ConnectedTopology ? max(2Nsubsteps + 3, old_halos[2]) : old_halos[2]

    new_halos = (Hx, Hy, old_halos[3])
    if new_halos == old_halos
        return grid
    else
        return with_halo(new_halos, grid)
    end
end

function materialize_solver(mom::SplitExplicitMomentumEquation, grid)
    new_auxiliaries = Auxiliaries(mom.rheology, grid)
    new_solver      = SplitExplicitSolver(grid; substeps = mom.solver.substeps)
    new_stress      = (bottom = materialize_stress(mom.external_momentum_stresses.bottom, grid),
                       top    = materialize_stress(mom.external_momentum_stresses.top, grid))

    return SeaIceMomentumEquation(mom.coriolis,
                                  mom.rheology,
                                  new_auxiliaries,
                                  new_solver,
                                  mom.free_drift,
                                  new_stress,
                                  mom.minimum_concentration,
                                  mom.minimum_mass)
end

# Reset the velocities to the previous time step
# This does nothing for a FE model, but is necessary for an RK model.
reset_velocities!(u, v, timestepper) = nothing

function reset_velocities!(u, v, timestepper::SplitRungeKuttaTimeStepper)
    parent(u) .= parent(timestepper.Ψ⁻.u)
    parent(v) .= parent(timestepper.Ψ⁻.v)
    return nothing
end

"""
    time_step_momentum!(model, dynamics::SplitExplicitMomentumEquation, Δt)

function for stepping u and v in the case of _explicit_ solvers.
The sea-ice momentum equations are characterized by smaller time-scale than
sea-ice ice_thermodynamics and sea-ice tracer advection, therefore explicit rheologies require
substepping over a set number of substeps.
"""
function time_step_momentum!(model, dynamics::SplitExplicitMomentumEquation, Δt)

    grid = model.velocities.u.grid
    arch = architecture(grid)

    # Unwrap variables
    rheology      = dynamics.rheology
    u, v          = model.velocities
    free_drift    = dynamics.free_drift
    clock         = model.clock
    coriolis      = dynamics.coriolis
    massmin       = dynamics.minimum_mass
    ℵmin          = dynamics.minimum_concentration
    u_forcing     = model.forcing.u
    v_forcing     = model.forcing.v
    Gu            = model.timestepper.Gⁿ.u
    Gv            = model.timestepper.Gⁿ.v
    u_immersed_bc = u.boundary_conditions.immersed
    v_immersed_bc = v.boundary_conditions.immersed
    top_stress    = dynamics.external_momentum_stresses.top
    bottom_stress = dynamics.external_momentum_stresses.bottom
    model_fields  = merge(dynamics.auxiliaries.fields, model.velocities,
                       (; h = model.ice_thickness,
                          ℵ = model.ice_concentration,
                          ρ = model.sea_ice_density))

    initialize_rheology!(model, dynamics.rheology)
    reset_velocities!(u, v, model.timestepper)

    params = dynamics.solver.kernel_parameters

    tendencies_kernel!, _ = configure_kernel(arch, grid, params, _compute_velocity_tendencies!)
    u_velocity_kernel!, _ = configure_kernel(arch, grid, params, _u_velocity_step!)
    v_velocity_kernel!, _ = configure_kernel(arch, grid, params, _v_velocity_step!)

    substeps = dynamics.solver.substeps

    u_args = (u, Gu, grid, Δt, substeps, rheology, model_fields, free_drift, clock, massmin, ℵmin, top_stress, bottom_stress)
    v_args = (v, Gv, grid, Δt, substeps, rheology, model_fields, free_drift, clock, massmin, ℵmin, top_stress, bottom_stress)
    G_args = (Gu, Gv, grid, Δt, rheology, model_fields, clock, coriolis, u_immersed_bc, v_immersed_bc, top_stress, bottom_stress, model.forcing)

    u_fill_halo_args = (u.data, u.boundary_conditions, u.indices, instantiated_location(u), grid, u.communication_buffers)
    v_fill_halo_args = (v.data, v.boundary_conditions, v.indices, instantiated_location(v), grid, v.communication_buffers)
    stresses_args    = (model_fields, grid, rheology, Δt)

    GC.@preserve G_args v_args u_args u_fill_halo_args v_fill_halo_args stresses_args begin
        # We need to timestep ~150 substeps, which means
        # launching ~1000 very small kernels: we are limited by
        # latency of argument conversion to GPU-compatible values.
        # To alleviate this penalty we convert first and then we substep!
        converted_u_args = convert_to_device(arch, u_args)
        converted_v_args = convert_to_device(arch, v_args)
        converted_G_args = convert_to_device(arch, G_args)

        # Do not convert args for fill halo regions if we are in a distributed scenario
        # (We need to know that we are passing a `DistributedGrid`)
        if arch isa Distributed
            converted_u_halo = u_fill_halo_args
            converted_v_halo = v_fill_halo_args
        else
            converted_u_halo = convert_to_device(arch, u_fill_halo_args)
            converted_v_halo = convert_to_device(arch, v_fill_halo_args)
        end

        converted_stresses_args = convert_to_device(arch, stresses_args)

        for substep in 1 : substeps
            fill_halo_regions!(converted_u_halo...; only_local_halos = true)
            fill_halo_regions!(converted_v_halo...; only_local_halos = true)

            # Compute stresses! depending on the particular rheology implementation
            compute_stresses!(dynamics, converted_stresses_args...)
            tendencies_kernel!(converted_G_args...)
            u_velocity_kernel!(converted_u_args...)
            v_velocity_kernel!(converted_v_args...)
        end
    end

    finalize_rheology!(model_fields, rheology)

    return nothing
end

@kernel function _u_velocity_step!(u, Gu, grid, Δt, substeps, rheology,
                                   fields, free_drift, clock,
                                   minimum_mass, minimum_concentration,
                                   u_top_stress, u_bottom_stress)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    mᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)
    ℵᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, fields.ℵ)

    Δτ = compute_substep_Δtᶠᶜᶜ(i, j, grid, Δt, rheology, substeps, fields)

    # Implicit part of the stress that depends linearly on the velocity
    τuᵢ = ( implicit_τx_coefficient(i, j, kᴺ, grid, u_bottom_stress, clock, fields)
          - implicit_τx_coefficient(i, j, kᴺ, grid, u_top_stress, clock, fields)) / mᵢ * ℵᵢ

    τuᵢ = ifelse(mᵢ ≤ 0, zero(grid), τuᵢ)
    uᴰ  = @inbounds (u[i, j, 1] + Δτ * Gu[i, j, 1]) / (1 + Δτ * τuᵢ) # dynamical velocity
    uᶠ  = free_drift_u(i, j, kᴺ, grid, free_drift, clock, fields) # free drift velocity

    # If the ice mass or the ice concentration are below a certain threshold,
    # the sea ice velocity is set to the free drift velocity. If no ice is
    # present above roundoff, the sea ice velocity is set to zero.
    marginal_ice = (mᵢ > eps(typeof(mᵢ))) & (ℵᵢ > eps(typeof(ℵᵢ)))
    active_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)
    active  = !peripheral_node(i, j, kᴺ, grid, Face(), Center(), Center())

    @inbounds u[i, j, 1] = ifelse(active_ice, uᴰ, ifelse(marginal_ice, uᶠ, zero(grid))) * active
end

@kernel function _v_velocity_step!(v, Gv, grid, Δt, substeps, rheology,
                                   fields, free_drift, clock,
                                   minimum_mass, minimum_concentration,
                                   v_top_stress, v_bottom_stress)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    mᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)
    ℵᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, fields.ℵ)

    Δτ = compute_substep_Δtᶜᶠᶜ(i, j, grid, Δt, rheology, substeps, fields)

    # Implicit part of the stress that depends linearly on the velocity
    τvᵢ = ( implicit_τy_coefficient(i, j, kᴺ, grid, v_bottom_stress, clock, fields)
          - implicit_τy_coefficient(i, j, kᴺ, grid, v_top_stress, clock, fields)) / mᵢ * ℵᵢ

    τvᵢ = ifelse(mᵢ ≤ 0, zero(grid), τvᵢ)

    vᴰ = @inbounds (v[i, j, 1] + Δτ * Gv[i, j, 1]) / (1 + Δτ * τvᵢ)# dynamical velocity
    vᶠ = free_drift_v(i, j, kᴺ, grid, free_drift, clock, fields)  # free drift velocity

    # If the ice mass or the ice concentration are below a certain threshold,
    # the sea ice velocity is set to the free drift velocity. If no ice is
    # present above roundoff, the sea ice velocity is set to zero.
    marginal_ice = (mᵢ > eps(typeof(mᵢ))) & (ℵᵢ > eps(typeof(ℵᵢ)))
    active_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)
    active  = !peripheral_node(i, j, kᴺ, grid, Center(), Face(), Center())

    @inbounds v[i, j, 1] = ifelse(active_ice, vᴰ, ifelse(marginal_ice, vᶠ, zero(grid))) * active
end
