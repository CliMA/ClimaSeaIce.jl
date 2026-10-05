using Oceananigans: instantiated_location
using Oceananigans.Architectures: convert_to_device, architecture, on_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.DistributedComputations: DistributedGrid
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, immersed_cell
using Oceananigans.Grids: AbstractGrid, halo_size, topology, with_halo, peripheral_node,
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

struct SplitExplicitSolver{I, K, A}
    substeps :: I
    kernel_parameters :: K
    active_cells :: A
end

SplitExplicitSolver(substeps, kernel_parameters) = SplitExplicitSolver(substeps, kernel_parameters, nothing)

"""
    SplitExplicitSolver(grid::AbstractGrid; substeps=120)

Creates a `SplitExplicitSolver` that controls the dynamical evolution of sea-ice momentum
by subcycling `substeps` times in between each ice_thermodynamics / tracer advection time step.

The default number of substeps is 120.

On an `ImmersedBoundaryGrid` with an active-columns map (`active_cells_map = true` or
`active_z_columns = true`), the substeps are computed only near columns that are not immersed;
see [`ImmersedActiveCells`](@ref).
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

maybe_extended_grid(mom::SplitExplicitMomentumEquation, grid::DistributedGrid) = maybe_extended_grid(mom.solver, grid)

# Halo-extended velocity grid for a split-explicit solver on a distributed grid.
function maybe_extended_grid(solver::SplitExplicitSolver, grid::DistributedGrid)
    old_halos = halo_size(grid)
    Nsubsteps = solver.substeps
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

#####
##### Skipping immersed columns
#####

"""
    ImmersedActiveCells

Index lists (and the kernels launched over them) that restrict the split-explicit substeps
to the points near columns that are not immersed. Built once, when the solver is materialized,
for an `ImmersedBoundaryGrid` with an active-columns map.

A point `(i, j)` is kept when any of the cells `(i-1:i, j-1:j)` is not immersed. These are the cells that set
- the `(Center, Center)` point `(i, j)`,
- the `(Face, Face)` point `(i, j)`, whose mass is interpolated from the four cells,
- the `(Face, Center)` point `(i, j)`, whose mass is interpolated from cells `i-1` and `i`, and
- the `(Center, Face)` point `(i, j)`, whose mass is interpolated from cells `j-1` and `j`.

Immersed columns hold no ice (`update_state!` masks the prognostic fields there), so at every point that is not kept
the kernels would only ever write the values the fields start with: zero velocities (the velocity
points are peripheral), unchanged (zero) stresses, zero viscosities and the relaxation parameter
`max_relaxation_parameter`. Skipping these points therefore gives bit-for-bit the same velocities,
stresses and relaxation parameter. The viscosities `ζ` and the diagnostic `Δ` of the
`ElastoViscoPlasticRheology` are not computed at the skipped points (where the full kernels can
produce NaNs in the halos); they are only ever read at the point where they were just computed.

The lists are sorted, so that consecutive entries are neighbours in memory, and hold signed
indices because the stress range extends into the halos.
"""
struct ImmersedActiveCells{L, K}
    lists :: L   # (; velocity, stress) index lists, on the architecture
    kernels :: K # (; u, v, stress) kernels launched over the lists
end

# Velocity kernels run over their `kernel_parameters`, stress kernels over the rheology-specific range
kernel_ranges(grid, ::Symbol) = (1:size(grid, 1), 1:size(grid, 2))
kernel_ranges(grid, ::KernelParameters{S, O}) where {S, O} = Tuple(1+o:s+o for (s, o) in zip(S, O))

# Only grids with an active-columns map skip immersed columns
immersed_active_cells(grid, rheology, auxiliaries, kernel_parameters) = nothing

function immersed_active_cells(grid::ImmersedBoundaryGrid, rheology, auxiliaries, kernel_parameters)
    isnothing(grid.active_z_columns) && return nothing

    arch = architecture(grid)
    velocity_list = active_cells_list(grid, kernel_ranges(grid, kernel_parameters))

    stress_ranges = stress_kernel_ranges(rheology, grid)
    stress_list = isnothing(stress_ranges) ? nothing : active_cells_list(grid, stress_ranges)

    u = configure_mapped_kernel(arch, grid, _u_velocity_step!, velocity_list)
    v = configure_mapped_kernel(arch, grid, _v_velocity_step!, velocity_list)
    stress = isnothing(stress_list) ? auxiliaries.kernels :
             mapped_stress_kernels(auxiliaries.kernels, rheology, arch, grid, stress_list)

    lists = (; velocity = velocity_list, stress = stress_list)
    return ImmersedActiveCells(lists, (; u, v, stress))
end

@kernel function _compute_active_neighbourhood!(mask, grid, i₀, j₀)
    i′, j′ = @index(Global, NTuple)
    i  = i′ + i₀
    j  = j′ + j₀
    kᴺ = size(grid, 3)

    active = !immersed_cell(i,   j,   kᴺ, grid) |
             !immersed_cell(i-1, j,   kᴺ, grid) |
             !immersed_cell(i,   j-1, kᴺ, grid) |
             !immersed_cell(i-1, j-1, kᴺ, grid)

    @inbounds mask[i′, j′] = active
end

# The `(i, j)` indices in `ranges` near at least one column that is not immersed
function active_cells_list(grid, ranges)
    arch = architecture(grid)
    i₀ = first(ranges[1]) - 1
    j₀ = first(ranges[2]) - 1

    mask = on_architecture(arch, zeros(Bool, length.(ranges)...))
    kernel!, _ = configure_kernel(arch, grid, KernelParameters(size(mask), (0, 0)), _compute_active_neighbourhood!)
    kernel!(mask, grid, i₀, j₀)

    # Built once, on the CPU; `findall` returns the indices sorted with `i` fastest
    indices = findall(on_architecture(CPU(), mask))
    list = [(Int32(I[1] + i₀), Int32(I[2] + j₀)) for I in indices]

    return on_architecture(arch, list)
end

function materialize_solver(mom::SplitExplicitMomentumEquation, grid)
    new_auxiliaries  = Auxiliaries(mom.rheology, grid)
    kernel_parameters = SplitExplicitSolver(grid; substeps = mom.solver.substeps).kernel_parameters
    active_cells     = immersed_active_cells(grid, mom.rheology, new_auxiliaries, kernel_parameters)
    new_solver       = SplitExplicitSolver(mom.solver.substeps, kernel_parameters, active_cells)
    new_basal_stress = materialize_basal_stress(mom.basal_stress, grid)
    new_free_surface = materialize_free_surface(mom.free_surface.η₀, mom.free_surface.g, grid)
    new_stress       = (bottom = materialize_stress(mom.external_momentum_stresses.bottom, grid),
                        top    = materialize_stress(mom.external_momentum_stresses.top, grid))

    # Repoint the free drift at the re-gridded stresses so it shares the same (extended) fields
    new_free_drift  = materialize_free_drift(mom.free_drift, new_stress.top, new_stress.bottom)

    return SeaIceMomentumEquation(mom.coriolis,
                                  mom.rheology,
                                  new_auxiliaries,
                                  new_solver,
                                  new_free_drift,
                                  new_stress,
                                  new_basal_stress,
                                  new_free_surface,
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
    basal_stress  = dynamics.basal_stress
    bottom_stress = dynamics.external_momentum_stresses.bottom
    model_fields  = merge(dynamics.auxiliaries.fields, model.velocities,
                       (; h = model.ice_thickness,
                          ℵ = model.ice_concentration,
                          ρ = model.sea_ice_density,
                          free_surface = dynamics.free_surface))

    reset_velocities!(u, v, model.timestepper)
    initialize_rheology!(model, dynamics.rheology)

    # Refresh the externally-provided stresses / velocities and fill their (extended) halos once per time step.
    update_external_stress!(top_stress, grid)
    update_external_stress!(bottom_stress, grid)
    update_free_surface!(dynamics.free_surface)

    params = dynamics.solver.kernel_parameters

    u_velocity_kernel!, _ = configure_kernel(arch, grid, params, _u_velocity_step!)
    v_velocity_kernel!, _ = configure_kernel(arch, grid, params, _v_velocity_step!)

    substeps = dynamics.solver.substeps

    u_args = (u, grid, Δt, substeps, rheology, model_fields, free_drift, clock, coriolis, massmin, ℵmin, u_immersed_bc, top_stress, bottom_stress, basal_stress, u_forcing)
    v_args = (v, grid, Δt, substeps, rheology, model_fields, free_drift, clock, coriolis, massmin, ℵmin, v_immersed_bc, top_stress, bottom_stress, basal_stress, v_forcing)

    u_fill_halo_args = (u.data, u.boundary_conditions, u.indices, instantiated_location(u), grid, u.communication_buffers)
    v_fill_halo_args = (v.data, v.boundary_conditions, v.indices, instantiated_location(v), grid, v.communication_buffers)
    stresses_args    = (model_fields, grid, rheology, Δt, u_immersed_bc, v_immersed_bc)

    GC.@preserve v_args u_args u_fill_halo_args v_fill_halo_args stresses_args begin
        # We need to timestep ~150 substeps, which means
        # launching ~1000 very small kernels: we are limited by
        # latency of argument conversion to GPU-compatible values.
        # To alleviate this penalty we convert first and then we substep!
        converted_u_args = convert_to_device(arch, u_args)
        converted_v_args = convert_to_device(arch, v_args)

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

        fill_halo_regions!(converted_u_halo...; only_local_halos = true)
        fill_halo_regions!(converted_v_halo...; only_local_halos = true)

        kernels = substep_kernels(dynamics.solver.active_cells, dynamics.auxiliaries.kernels,
                                  u_velocity_kernel!, v_velocity_kernel!)

        for substep in 1 : substeps
            momentum_substep!(kernels, substep, converted_stresses_args, converted_u_args, converted_v_args,
                              converted_u_halo, converted_v_halo)
        end
    end

    finalize_rheology!(model_fields, rheology)

    return nothing
end

# Kernels over the whole domain, or only near the columns that are not immersed
substep_kernels(::Nothing, stress_kernels, u_kernel!, v_kernel!) = (; stress = stress_kernels, u = u_kernel!, v = v_kernel!)
substep_kernels(active_cells::ImmersedActiveCells, args...) = active_cells.kernels

function momentum_substep!(kernels, substep, stresses_args, u_args, v_args, u_halo, v_halo)
    # Compute stresses! depending on the particular rheology implementation
    compute_stresses!(kernels.stress, stresses_args...)

    # Alternating leap-frog.
    if iseven(substep)
        kernels.u(u_args...)
        fill_halo_regions!(u_halo...; only_local_halos = true)
        kernels.v(v_args...)
        fill_halo_regions!(v_halo...; only_local_halos = true)
    else
        kernels.v(v_args...)
        fill_halo_regions!(v_halo...; only_local_halos = true)
        kernels.u(u_args...)
        fill_halo_regions!(u_halo...; only_local_halos = true)
    end

    return nothing
end

@kernel function _u_velocity_step!(u, grid, Δt, substeps, rheology,
                                   fields, free_drift, clock, coriolis,
                                   minimum_mass, minimum_concentration,
                                   u_immersed_bc, u_top_stress, u_bottom_stress, basal_stress, u_forcing)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    mᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)
    ℵᵢ = ℑxᶠᵃᵃ(i, j, kᴺ, grid, fields.ℵ)

    Δτ = compute_substep_Δtᶠᶜᶜ(i, j, grid, Δt, rheology, substeps, fields)

    Gu = u_velocity_tendency(i, j, grid, Δτ, rheology, fields, clock, coriolis,
                             u_immersed_bc, u_top_stress, u_bottom_stress, u_forcing)

    # The external stresses act over the ice-covered fraction; the basal stress carries its own
    τuᵢ = ( implicit_τx_coefficient(i, j, kᴺ, grid, u_bottom_stress, clock, fields) / mᵢ * ℵᵢ
          - implicit_τx_coefficient(i, j, kᴺ, grid, u_top_stress, clock, fields) / mᵢ * ℵᵢ
          + basal_τx_coefficient(i, j, kᴺ, grid, basal_stress, fields) / mᵢ)

    τuᵢ = ifelse(mᵢ ≤ 0, zero(grid), τuᵢ)
    uᴰ  = @inbounds (u[i, j, 1] + Δτ * Gu) / (1 + Δτ * τuᵢ) # dynamical velocity
    uᶠ  = free_drift_u(i, j, kᴺ, grid, free_drift, clock, fields) # free drift velocity

    # If the ice mass or the ice concentration are below a certain threshold,
    # the sea ice velocity is set to the free drift velocity. If no ice is
    # present above roundoff, the sea ice velocity is set to zero.
    marginal_ice = (mᵢ > eps(typeof(mᵢ))) & (ℵᵢ > eps(typeof(ℵᵢ)))
    active_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)
    active  = !peripheral_node(i, j, kᴺ, grid, Face(), Center(), Center())

    @inbounds u[i, j, 1] = ifelse(active_ice, uᴰ, ifelse(marginal_ice, uᶠ, zero(grid))) * active
end

@kernel function _v_velocity_step!(v, grid, Δt, substeps, rheology,
                                   fields, free_drift, clock, coriolis,
                                   minimum_mass, minimum_concentration,
                                   v_immersed_bc, v_top_stress, v_bottom_stress, basal_stress, v_forcing)

    i, j = @index(Global, NTuple)
    kᴺ   = size(grid, 3)

    mᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, ice_mass, fields.h, fields.ℵ, fields.ρ)
    ℵᵢ = ℑyᵃᶠᵃ(i, j, kᴺ, grid, fields.ℵ)

    Δτ = compute_substep_Δtᶜᶠᶜ(i, j, grid, Δt, rheology, substeps, fields)

    Gv = v_velocity_tendency(i, j, grid, Δτ, rheology, fields, clock, coriolis,
                             v_immersed_bc, v_top_stress, v_bottom_stress, v_forcing)

    # Implicit part of the stress that depends linearly on the velocity
    τvᵢ = ( implicit_τy_coefficient(i, j, kᴺ, grid, v_bottom_stress, clock, fields) / mᵢ * ℵᵢ
          - implicit_τy_coefficient(i, j, kᴺ, grid, v_top_stress, clock, fields) / mᵢ * ℵᵢ
          + basal_τy_coefficient(i, j, kᴺ, grid, basal_stress, fields) / mᵢ)

    τvᵢ = ifelse(mᵢ ≤ 0, zero(grid), τvᵢ)

    vᴰ = @inbounds (v[i, j, 1] + Δτ * Gv) / (1 + Δτ * τvᵢ)# dynamical velocity
    vᶠ = free_drift_v(i, j, kᴺ, grid, free_drift, clock, fields)  # free drift velocity

    # If the ice mass or the ice concentration are below a certain threshold,
    # the sea ice velocity is set to the free drift velocity. If no ice is
    # present above roundoff, the sea ice velocity is set to zero.
    marginal_ice = (mᵢ > eps(typeof(mᵢ))) & (ℵᵢ > eps(typeof(ℵᵢ)))
    active_ice = (mᵢ ≥ minimum_mass) & (ℵᵢ ≥ minimum_concentration)
    active  = !peripheral_node(i, j, kᴺ, grid, Center(), Face(), Center())

    @inbounds v[i, j, 1] = ifelse(active_ice, vᴰ, ifelse(marginal_ice, vᶠ, zero(grid))) * active
end
