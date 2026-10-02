using Oceananigans: instantiated_location
using Oceananigans.Architectures: convert_to_device, architecture, device, on_architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.DistributedComputations: DistributedGrid
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

struct SplitExplicitSolver{I, K, S}
    substeps :: I
    kernel_parameters :: K
    ice_free_cells :: S
end

SplitExplicitSolver(substeps, kernel_parameters) = SplitExplicitSolver(substeps, kernel_parameters, nothing)

"""
    SplitExplicitSolver(grid::AbstractGrid; substeps=120, skip_ice_free_cells=false)

Creates a `SplitExplicitSolver` that controls the dynamical evolution of sea-ice momentum
by subcycling `substeps` times in between each ice_thermodynamics / tracer advection time step.

The default number of substeps is 120.

If `skip_ice_free_cells = true`, all substeps after the first are computed only in cells
that have ice in them or in a neighbouring cell. The answer does not change, since
the ice-free cells already stop changing after the first substep. See [`SkipIceFreeCells`](@ref).
"""
SplitExplicitSolver(grid::AbstractGrid; substeps=120, skip_ice_free_cells=false) =
    SplitExplicitSolver(substeps, :xy, ice_free_cells_option(skip_ice_free_cells))

# When no grid is provided, we assume a serial grid with default kernel parameters
SplitExplicitSolver(; substeps=120, skip_ice_free_cells=false) =
    SplitExplicitSolver(substeps, :xy, ice_free_cells_option(skip_ice_free_cells))

const SplitExplicitMomentumEquation = SeaIceMomentumEquation{<:SplitExplicitSolver}

# Shenanigans for extending the halos in Distributed grids

function SplitExplicitSolver(grid::DistributedGrid; substeps=120, skip_ice_free_cells=false)
    Nx, Ny, _ = size(grid)
    Hx, Hy, _ = halo_size(grid)
    TX, TY, _ = topology(grid)
    kernel_sizes = map(split_explicit_kernel_size, (TX, TY), (Nx, Ny), (Hx, Hy))
    return SplitExplicitSolver(substeps, KernelParameters(kernel_sizes...), ice_free_cells_option(skip_ice_free_cells))
end

skips_ice_free_cells(solver::SplitExplicitSolver) = !isnothing(solver.ice_free_cells)

Base.summary(solver::SplitExplicitSolver) =
    string("SplitExplicitSolver(substeps=", solver.substeps,
           ", skip_ice_free_cells=", skips_ice_free_cells(solver), ")")

#####
##### Skipping ice-free cells
#####

"""
    SkipIceFreeCells

Lets the `SplitExplicitSolver` skip the points where nothing changes during the substeps.

A point `(i, j)` is "iced" when any of the cells `(i-1:i, j-1:j)` has nonzero ice mass. These are the cells that set
- the `(Center, Center)` point `(i, j)`,
- the `(Face, Face)` point `(i, j)`, whose mass is interpolated from the four cells,
- the `(Face, Center)` point `(i, j)`, whose mass is interpolated from cells `i-1` and `i`, and
- the `(Center, Face)` point `(i, j)`, whose mass is interpolated from cells `j-1` and `j`.

At points that are not iced, every substep writes the same values: zero velocities, unchanged
stresses (zero mass), zero viscosities (zero ice strength) and the relaxation parameter of an
ice-free cell. A point that is not iced and already holds these values is "settled", and launching
the kernels there changes nothing.

At the beginning of each `time_step_momentum!`, two lists of `(i, j)` indices are built for each kernel range:
- the points that are iced or not yet settled (for example, where the ice has just melted), which the first substep computes, and
- the iced points, which all the other substeps compute: after the first substep every other point is settled.

The ice thickness and concentration do not change during the substeps, so neither do the lists.
The velocities and stresses are bit-for-bit the same as without skipping. The only exception is the
diagnostic `Δ` of the `ElastoViscoPlasticRheology`, which is left stale where there is no ice.
It is only ever read at the point where it was just computed.

The lists are filled in place, in buffers allocated once, so no memory is allocated during the time step.
The lists are built row by row: each row `j` is first counted, the row counts are summed to find where
each row starts in the list, and then each row is written in order of `i`. The indices are therefore
sorted, so that consecutive entries are neighbours in memory.
"""
struct SkipIceFreeCells{V, S}
    velocity :: V # `IceCoverMaps` over the velocity kernel range
    stress :: S   # `IceCoverMaps` over the stress kernel range (`nothing` if the rheology has no stress kernels)
end

# Before the solver knows its grid
SkipIceFreeCells() = SkipIceFreeCells(nothing, nothing)

ice_free_cells_option(skip::Bool) = skip ? SkipIceFreeCells() : nothing

struct IceCoverMaps{I, C, T, H, R, K}
    first_substep :: I # Indices computed during the first substep
    iced :: I          # Indices computed during all the other substeps
    row_counts :: C    # Number of indices of each list in each row
    row_offsets :: C   # Where each row starts in each list
    counts :: T        # Number of indices in each list, on the architecture
    host_counts :: H   # Number of indices in each list, on the CPU
    ranges :: R        # The `(i, j)` ranges covered by the lists
    kernels :: K       # Kernels counting the rows, summing the counts, and filling the rows
end

function IceCoverMaps(grid, ranges)
    arch = architecture(grid)
    dev  = device(arch)
    N    = prod(length.(ranges))
    Ny   = length(ranges[2])

    # The indices are signed because the stress range extends into the halos
    first_substep = on_architecture(arch, Vector{NTuple{2, Int32}}(undef, N))
    iced          = on_architecture(arch, Vector{NTuple{2, Int32}}(undef, N))
    row_counts    = on_architecture(arch, zeros(Int32, 2, Ny))
    row_offsets   = on_architecture(arch, zeros(Int32, 2, Ny))
    counts        = on_architecture(arch, zeros(Int32, 2))
    host_counts   = zeros(Int32, 2)

    workgroup = min(Ny, 64)
    kernels = (count = _count_ice_cover_rows!(dev, workgroup, Ny),
               sum   = _sum_ice_cover_rows!(dev, 1, 1), # The number of rows is small: a single work item
               fill  = _fill_ice_cover_rows!(dev, workgroup, Ny))

    return IceCoverMaps(first_substep, iced, row_counts, row_offsets, counts, host_counts, ranges, kernels)
end

# Velocity kernels run over their `kernel_parameters`, stress kernels over the rheology-specific range
kernel_ranges(grid, ::Symbol) = (1:size(grid, 1), 1:size(grid, 2))
kernel_ranges(grid, ::KernelParameters{S, O}) where {S, O} = Tuple(1+o:s+o for (s, o) in zip(S, O))

materialize_ice_free_cells(::Nothing, grid, rheology, kernel_parameters) = nothing

function materialize_ice_free_cells(::SkipIceFreeCells, grid, rheology, kernel_parameters)
    velocity = IceCoverMaps(grid, kernel_ranges(grid, kernel_parameters))
    stress_ranges = stress_kernel_ranges(rheology, grid)
    stress = isnothing(stress_ranges) ? nothing : IceCoverMaps(grid, stress_ranges)
    return SkipIceFreeCells(velocity, stress)
end

# Whether `(i, j)` belongs to the first-substep list and to the iced list
@inline function ice_cover_lists(i, j, grid, rheology, fields)
    h = fields.h
    ℵ = fields.ℵ
    ρ = fields.ρ

    # `≠ 0` rather than `> 0` so that cells with negative or NaN mass are still computed
    has_ice = (ice_mass(i,   j,   1, grid, h, ℵ, ρ) != 0) |
              (ice_mass(i-1, j,   1, grid, h, ℵ, ρ) != 0) |
              (ice_mass(i,   j-1, 1, grid, h, ℵ, ρ) != 0) |
              (ice_mass(i-1, j-1, 1, grid, h, ℵ, ρ) != 0)

    moving = @inbounds !iszero(fields.u[i, j, 1]) | !iszero(fields.v[i, j, 1])
    unsettled = moving | unsettled_stresses(i, j, grid, rheology, fields)

    return has_ice | unsettled, has_ice
end

@kernel function _count_ice_cover_rows!(row_counts, grid, rheology, fields, irange, j₀)
    j′ = @index(Global, Linear)
    j  = j′ + j₀

    n₁ = Int32(0)
    n₂ = Int32(0)
    for i in irange
        first_substep, iced = ice_cover_lists(i, j, grid, rheology, fields)
        n₁ += first_substep
        n₂ += iced
    end

    @inbounds row_counts[1, j′] = n₁
    @inbounds row_counts[2, j′] = n₂
end

@kernel function _sum_ice_cover_rows!(row_offsets, counts, row_counts)
    n₁ = Int32(0)
    n₂ = Int32(0)
    @inbounds for j′ in 1:size(row_counts, 2)
        row_offsets[1, j′] = n₁
        row_offsets[2, j′] = n₂
        n₁ += row_counts[1, j′]
        n₂ += row_counts[2, j′]
    end

    @inbounds counts[1] = n₁
    @inbounds counts[2] = n₂
end

@kernel function _fill_ice_cover_rows!(first_substep_list, iced_list, row_offsets, grid, rheology, fields, irange, j₀)
    j′ = @index(Global, Linear)
    j  = j′ + j₀

    n₁ = @inbounds row_offsets[1, j′]
    n₂ = @inbounds row_offsets[2, j′]
    for i in irange
        first_substep, iced = ice_cover_lists(i, j, grid, rheology, fields)
        index = (Int32(i), Int32(j))
        if first_substep
            n₁ += Int32(1)
            @inbounds first_substep_list[n₁] = index
        end
        if iced
            n₂ += Int32(1)
            @inbounds iced_list[n₂] = index
        end
    end
end

# Fill the lists and return views of their filled parts.
# This requires copying the two counts to the CPU, but allocates no memory.
function update_ice_cover_maps!(maps::IceCoverMaps, grid, rheology, fields)
    irange, jrange = maps.ranges
    j₀ = first(jrange) - 1

    maps.kernels.count(maps.row_counts, grid, rheology, fields, irange, j₀)
    maps.kernels.sum(maps.row_offsets, maps.counts, maps.row_counts)
    maps.kernels.fill(maps.first_substep, maps.iced, maps.row_offsets, grid, rheology, fields, irange, j₀)
    copyto!(maps.host_counts, maps.counts)

    first_substep = view(maps.first_substep, 1:Int(maps.host_counts[1]))
    iced          = view(maps.iced,          1:Int(maps.host_counts[2]))

    return first_substep, iced
end

function restricted_kernels(full_kernels, arch, grid, rheology, velocity_map, stress_map)
    u_kernel! = configure_mapped_kernel(arch, grid, _u_velocity_step!, velocity_map)
    v_kernel! = configure_mapped_kernel(arch, grid, _v_velocity_step!, velocity_map)
    stress_kernels = isnothing(stress_map) ? full_kernels.stress :
                     mapped_stress_kernels(full_kernels.stress, rheology, arch, grid, stress_map)
    return (; stress = stress_kernels, u = u_kernel!, v = v_kernel!)
end

# Kernels for the first substep and for all the other substeps
substep_kernels(::Nothing, full_kernels, args...) = (full_kernels, full_kernels)

function substep_kernels(skip::SkipIceFreeCells, full_kernels, arch, grid, rheology, fields)
    velocity_first, velocity_iced = update_ice_cover_maps!(skip.velocity, grid, rheology, fields)

    if isnothing(skip.stress)
        stress_first = stress_iced = nothing
    else
        stress_first, stress_iced = update_ice_cover_maps!(skip.stress, grid, rheology, fields)
    end

    first_kernels = restricted_kernels(full_kernels, arch, grid, rheology, velocity_first, stress_first)
    later_kernels = restricted_kernels(full_kernels, arch, grid, rheology, velocity_iced,  stress_iced)

    return first_kernels, later_kernels
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

function materialize_split_explicit_solver(solver, rheology, grid)
    kernel_parameters = SplitExplicitSolver(grid; substeps = solver.substeps).kernel_parameters
    ice_free_cells = materialize_ice_free_cells(solver.ice_free_cells, grid, rheology, kernel_parameters)
    return SplitExplicitSolver(solver.substeps, kernel_parameters, ice_free_cells)
end

function materialize_solver(mom::SplitExplicitMomentumEquation, grid)
    new_auxiliaries  = Auxiliaries(mom.rheology, grid)
    new_solver       = materialize_split_explicit_solver(mom.solver, mom.rheology, grid)
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

        # Possibly restrict the kernels to the points that change during the substeps (see `SkipIceFreeCells`)
        full_kernels = (; stress = dynamics.auxiliaries.kernels, u = u_velocity_kernel!, v = v_velocity_kernel!)
        first_kernels, later_kernels = substep_kernels(dynamics.solver.ice_free_cells, full_kernels,
                                                       arch, grid, rheology, model_fields)

        for substep in 1 : substeps
            kernels = substep == 1 ? first_kernels : later_kernels
            momentum_substep!(kernels, substep, converted_stresses_args, converted_u_args, converted_v_args,
                              converted_u_halo, converted_v_halo)
        end
    end

    finalize_rheology!(model_fields, rheology)

    return nothing
end

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
