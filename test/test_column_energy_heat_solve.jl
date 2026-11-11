using ClimaSeaIce.SeaIceThermodynamics:
    BulkSalinityDiffusion,
    SeaIceColumnDiscretization,
    ConductiveTemperatureTransport,
    ExponentialShortwaveAbsorption,
    FixedDrainedIceSalinityProfile,
    FixedSalinityBrinePocketEnergyRelation,
    FluxBoundary,
    MaykutUntersteinerConductivity,
    NoSalinityTransport,
    PrescribedTemperature,
    QuadraticLiquidusEnergyRelation,
    assemble_column_energy_system!,
    column_energy_time_step!,
    column_salinity_time_step!,
    compute_column_shortwave_flux!,
    compute_column_thermodynamic_diagnostics!,
    compute_column_transport_coefficients!,
    evolving_salinity_mushy_thermodynamics,
    face_thermal_conductivity,
    ice_thermal_conductivity,
    internal_energy,
    prescribed_salinity_enthalpy_thermodynamics,
    solve_column_energy_system!,
    temperature
using Oceananigans
using Oceananigans.TimeSteppers: Clock
using Oceananigans.Fields: interior, set!, ConstantField
using Test

# Reference conservative remap of piecewise-constant layer averages; uncovered target intervals take `fill_value`.
function conservative_column_remap(source_values, source_faces, target_faces; fill_value = zero(eltype(source_values)))
    return map(1:length(target_faces)-1) do k
        zᵇ, zᵗ = target_faces[k], target_faces[k+1]
        overlaps = [max(min(zᵗ, source_faces[s+1]) - max(zᵇ, source_faces[s]), 0) for s in eachindex(source_values)]
        (fill_value * (zᵗ - zᵇ) + sum((source_values .- fill_value) .* overlaps)) / (zᵗ - zᵇ)
    end
end

# Set the SeaIceColumnDiscretization interfaces to a previous→current thickness change. `move = :top` keeps the
# base at `surface` and moves the surface (the default, e.g. surface melt/snow-ice); `move = :base` keeps the
# surface at `surface` and moves the base (basal congelation/melt).
function set_column_metric!(grid, previous, current; surface = 0, move = :top)
    z = grid.z
    if move === :top
        fill!(z.hb⁻, surface); fill!(z.hbⁿ, surface)
        fill!(z.hs⁻, surface + previous); fill!(z.hsⁿ, surface + current)
    else
        fill!(z.hs⁻, surface); fill!(z.hsⁿ, surface)
        fill!(z.hb⁻, surface - previous); fill!(z.hbⁿ, surface - current)
    end
    return nothing
end

@testset "CICE-compatible thermal conductivity on a column grid" begin
    maykut = MaykutUntersteinerConductivity(Float64)

    grid = RectilinearGrid(size = 4,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    profile = FixedDrainedIceSalinityProfile(Float64)
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    energy_transport = ConductiveTemperatureTransport(conductivity = maykut)
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = profile,
        energy_transport,
        heat_boundary_conditions = insulating_boundaries())

    set!(thermodynamics;
         bulk_salinity = profile,
         temperature = z -> -15 + 10z)

    compute_column_transport_coefficients!(thermodynamics)

    temperatures = column_values(thermodynamics.fields.temperature)
    salinities = column_values(thermodynamics.fields.bulk_salinity)
    face_conductivities = column_values(thermodynamics.auxiliary.thermal_conductivity)
    dz = 0.25

    @test isapprox(face_conductivities[1],
                   ice_thermal_conductivity(maykut, temperatures[1], salinities[1]);
                   atol = 1e-14, rtol = 0)
    @test isapprox(face_conductivities[end],
                   ice_thermal_conductivity(maykut, temperatures[end], salinities[end]);
                   atol = 1e-14, rtol = 0)

    for k in 2:4
        expected_K = face_thermal_conductivity(maykut,
                                               temperatures[k-1],
                                               salinities[k-1],
                                               dz,
                                               temperatures[k],
                                               salinities[k],
                                               dz)
        @test isapprox(face_conductivities[k], expected_K; atol = 1e-14, rtol = 0)
    end
end

function solve_reference_vector_tridiagonal(lower_diagonal, diagonal, upper_diagonal, rhs)
    n = length(rhs)
    solution = similar(rhs)
    scratch = zeros(eltype(rhs), max(n - 1, 0))
    β = first(diagonal)

    solution[1] = rhs[1] / β

    for k in 2:n
        scratch[k-1] = upper_diagonal[k-1] / β
        β = diagonal[k] - lower_diagonal[k-1] * scratch[k-1]
        solution[k] = (rhs[k] - lower_diagonal[k-1] * solution[k-1]) / β
    end

    for k in n-1:-1:1
        solution[k] -= scratch[k] * solution[k+1]
    end

    return solution
end

function pure_ice_conductive_column(N; diffusivity = 1.0)
    grid = RectilinearGrid(size = N,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_capacity = relation.phase_transitions.density *
                    relation.phase_transitions.heat_capacity
    conductivity = diffusivity * heat_capacity
    energy_transport = ConductiveTemperatureTransport(conductivity = conductivity)
    heat_boundary_conditions = insulating_boundaries()

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        heat_boundary_conditions)

    return thermodynamics, heat_capacity
end

function initialize_cosine_energy_mode!(thermodynamics, heat_capacity;
                                        base_energy = -2e7,
                                        amplitude = 1e5)
    set!(thermodynamics;
         bulk_salinity = 0.0,
         temperature = z -> (base_energy + amplitude * cos(pi * z)) / heat_capacity)

    return nothing
end

function spatial_manufactured_solution_error(N; dt = 0.01, diffusivity = 1.0)
    thermodynamics, heat_capacity = pure_ice_conductive_column(N; diffusivity)
    base_energy = -2e7
    amplitude = 1e5

    initialize_cosine_energy_mode!(thermodynamics, heat_capacity;
                                   base_energy,
                                   amplitude)

    column_step!(thermodynamics, insulating_fluxes(), dt)

    numerical_energy = column_values(thermodynamics.fields.internal_energy)
    exact_amplitude = amplitude / (1 + diffusivity * pi^2 * dt)
    expected_energy = [base_energy + exact_amplitude * cos(pi * (k - 0.5) / N)
                       for k in 1:N]

    return maximum(abs.(numerical_energy .- expected_energy))
end

function temporal_manufactured_solution_error(dt; N = 128,
                                              final_time = 0.16,
                                              diffusivity = 1.0)
    thermodynamics, heat_capacity = pure_ice_conductive_column(N; diffusivity)
    base_energy = -2e7
    amplitude = 1e5
    steps = round(Int, final_time / dt)

    initialize_cosine_energy_mode!(thermodynamics, heat_capacity;
                                   base_energy,
                                   amplitude)

    for _ in 1:steps
        column_step!(thermodynamics, insulating_fluxes(), dt)
    end

    numerical_energy = column_values(thermodynamics.fields.internal_energy)
    discrete_eigenvalue = 4 * N^2 * sin(pi / (2N))^2
    exact_amplitude = amplitude * exp(-diffusivity * discrete_eigenvalue * final_time)
    expected_energy = [base_energy + exact_amplitude * cos(pi * (k - 0.5) / N)
                       for k in 1:N]

    return maximum(abs.(numerical_energy .- expected_energy))
end

@testset "Column energy fixed-grid time step" begin
    grid = RectilinearGrid(size = 8,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_boundary_conditions = insulating_boundaries()
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0)
    dt = 5e3

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        heat_boundary_conditions)

    set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> -12 + 4z)

    initial_energy = column_integrated_energy(thermodynamics)
    initial_profile = column_values(thermodynamics.fields.internal_energy)

    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_transport_coefficients!(thermodynamics)
    assemble_column_energy_system!(thermodynamics, insulating_fluxes(), TEST_CLOCK, NO_FIELDS, CONSOLIDATED, dt)

    matrix = dense_column_energy_matrix(thermodynamics)
    rhs = column_values(thermodynamics.auxiliary.energy_rhs)
    expected_energy = matrix \ rhs

    solve_column_energy_system!(thermodynamics)
    stepped_energy = column_values(thermodynamics.fields.internal_energy)

    @test maximum(abs.(stepped_energy .- expected_energy)) < 1e-7
    @test abs(column_integrated_energy(thermodynamics) - initial_energy) < 1e-7
    @test maximum(abs.(stepped_energy .- initial_profile)) > 1e5

    constant = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        heat_boundary_conditions)

    set!(constant; bulk_salinity = 0.0, temperature = -10.0)
    constant_energy = column_integrated_energy(constant)
    constant_profile = column_values(constant.fields.internal_energy)

    column_step!(constant, insulating_fluxes(), dt)

    @test maximum(abs.(column_values(constant.fields.internal_energy) .- constant_profile)) < 1e-7
    @test abs(column_integrated_energy(constant) - constant_energy) < 1e-7
    @test maximum(abs.(column_values(constant.fields.temperature) .+ 10)) < 1e-12

    flux_fluxes = (top = 2.0, bottom = 0.5)
    flux_forced = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary()))

    set!(flux_forced; bulk_salinity = 0.0, temperature = -10.0)
    flux_forced_energy = column_integrated_energy(flux_forced)
    column_step!(flux_forced, flux_fluxes, 10.0)

    # Upward-positive fluxes: net energy gain is Δt·(Qᵇ − Qᵘ) = 10·(0.5 − 2.0).
    @test abs(column_integrated_energy(flux_forced) -
              flux_forced_energy -
              10.0 * (0.5 - 2.0)) < 1e-8

    shortwave = ExponentialShortwaveAbsorption(surface_transmission = 3.0,
                                               attenuation_scale = 0.25)
    shortwave_forced = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        shortwave_absorption = shortwave,
        heat_boundary_conditions)

    set!(shortwave_forced; bulk_salinity = 0.0, temperature = -10.0)
    shortwave_forced_energy = column_integrated_energy(shortwave_forced)
    compute_column_shortwave_flux!(shortwave_forced)

    @test column_values(shortwave_forced.auxiliary.shortwave_flux)[end] ≈ 3.0
    @test column_values(shortwave_forced.auxiliary.shortwave_flux)[1] ≈ 3.0 * exp(-4)

    column_step!(shortwave_forced, insulating_fluxes(), 100.0)
    expected_shortwave_change = 100.0 * 3.0 * (1 - exp(-4))

    @test abs(column_integrated_energy(shortwave_forced) - shortwave_forced_energy - expected_shortwave_change) < 1e-8

    bottom_temperature = -1.8
    top_flux = -3.0
    prescribed_bottom = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        heat_boundary_conditions = (top = FluxBoundary(),
                                    bottom = PrescribedTemperature(bottom_temperature)))

    set!(prescribed_bottom; bulk_salinity = 0.0, temperature = z -> -12 + 4z)

    initial_temperature = column_values(prescribed_bottom.fields.temperature)
    n = length(initial_temperature)
    Δz = 1 / n
    ρc = relation.phase_transitions.density * relation.phase_transitions.heat_capacity
    η = dt / (ρc * Δz)
    conductance = 2 * 2.0 / Δz
    interface_conductance = 2.0 / Δz
    lower_diagonal = zeros(n - 1)
    diagonal = ones(n)
    upper_diagonal = zeros(n - 1)
    rhs = copy(initial_temperature)

    for k in 1:n
        if k == 1
            diagonal[k] += η * (conductance + interface_conductance)
            upper_diagonal[k] = -η * interface_conductance
            rhs[k] += η * conductance * bottom_temperature
        elseif k == n
            lower_diagonal[k-1] = -η * interface_conductance
            diagonal[k] += η * interface_conductance
            rhs[k] -= η * top_flux  # upward-positive top flux enters the top cell as -Qᵘ
        else
            lower_diagonal[k-1] = -η * interface_conductance
            upper_diagonal[k] = -η * interface_conductance
            diagonal[k] += 2 * η * interface_conductance
        end
    end

    expected_temperature =
        solve_reference_vector_tridiagonal(lower_diagonal, diagonal, upper_diagonal, rhs)

    column_step!(prescribed_bottom, (top = top_flux, bottom = 0.0), dt)

    @test maximum(abs.(column_values(prescribed_bottom.fields.temperature) .-
                       expected_temperature)) < 1e-11
end

@testset "Nonlinear column conduction at one 90-minute ice step" begin
    Nz = 8
    Δz = 1 / Nz
    Δt = 90 * 60.0
    grid = RectilinearGrid(Float64; size = Nz, z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation, salinity_profile = FixedDrainedIceSalinityProfile(Float64),
        energy_transport = ConductiveTemperatureTransport(conductivity = MaykutUntersteinerConductivity(Float64)),
        heat_boundary_conditions = insulating_boundaries())
    set!(thermodynamics; temperature = z -> -2 - 18z)

    initial_energy = column_values(thermodynamics.fields.internal_energy)
    initial_temperature = column_values(thermodynamics.fields.temperature)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_transport_coefficients!(thermodynamics)
    assemble_column_energy_system!(thermodynamics, insulating_fluxes(), TEST_CLOCK, NO_FIELDS, CONSOLIDATED, Δt)

    # Each side of a conductive face has its own dT/dE. A face-averaged value
    # does not give the Jacobian of the two-layer heat flux near melting.
    matrix = dense_column_energy_matrix(thermodynamics)
    K = column_values(thermodynamics.auxiliary.thermal_conductivity)
    TE = column_values(thermodynamics.fields.temperature_energy_derivative)
    for k in 1:Nz-1
        @test isapprox(matrix[k, k+1], -Δt * K[k+1] * TE[k+1] / Δz^2; rtol = 1e-12)
        @test isapprox(matrix[k+1, k], -Δt * K[k+1] * TE[k] / Δz^2; rtol = 1e-12)
    end

    column_step!(thermodynamics, insulating_fluxes(), Δt)
    final_energy = column_values(thermodynamics.fields.internal_energy)
    final_temperature = column_values(thermodynamics.fields.temperature)
    compute_column_transport_coefficients!(thermodynamics)
    K = column_values(thermodynamics.auxiliary.thermal_conductivity)
    F = zeros(Nz + 1)
    for k in 2:Nz
        F[k] = K[k] * (final_temperature[k] - final_temperature[k-1]) / Δz
    end
    residual = [final_energy[k] - initial_energy[k] - Δt * (F[k+1] - F[k]) / Δz for k in 1:Nz]
    @test maximum(abs, residual) < 20 # J/m³, versus O(10⁶) J/m³ conductive tendencies
    @test minimum(final_temperature) >= minimum(initial_temperature) - 1e-5
    @test maximum(final_temperature) <= maximum(initial_temperature) + 1e-5
    @test isapprox(sum(final_energy), sum(initial_energy); atol = 1e-3)
end

@testset "Column energy moving-grid metric" begin
    grid = RectilinearGrid(size = 4,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_boundary_conditions = insulating_boundaries()
    no_transport = ConductiveTemperatureTransport(conductivity = 0.0)

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = no_transport,
        heat_boundary_conditions)

    set!(thermodynamics; bulk_salinity = 0.0, temperature = -10.0)

    initial_energy = column_integrated_energy(thermodynamics)
    initial_profile = column_values(thermodynamics.fields.internal_energy)
    current_metric = 0.75

    set_column_metric!(grid, 1.0, current_metric)
    column_step!(thermodynamics, insulating_fluxes(), 10.0)

    @test column_integrated_energy(thermodynamics) ≈ current_metric * initial_energy
    @test column_values(thermodynamics.fields.internal_energy) ≈ initial_profile

    set_column_metric!(grid, 1.0, 1.0)
    set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> -12 + 8z)
    nonuniform_initial_profile = column_values(thermodynamics.fields.internal_energy)
    source_faces = collect(range(0.0, 1.0; length = 5))
    contracted_faces = collect(range(0.0, 0.75; length = 5))
    expected_contracted_profile =
        conservative_column_remap(nonuniform_initial_profile, source_faces, contracted_faces)

    set_column_metric!(grid, 1.0, 0.75)
    column_step!(thermodynamics, insulating_fluxes(), 10.0)

    @test column_values(thermodynamics.fields.internal_energy) ≈ expected_contracted_profile

    set_column_metric!(grid, 1.0, 1.0)
    set!(thermodynamics; bulk_salinity = 0.0, temperature = z -> -12 + 8z)
    nonuniform_initial_profile = column_values(thermodynamics.fields.internal_energy)
    expanded_faces = collect(range(0.0, 1.25; length = 5))
    expected_expanded_profile =
        conservative_column_remap(nonuniform_initial_profile,
                                  source_faces,
                                  expanded_faces;
                                  fill_value = last(nonuniform_initial_profile))

    set_column_metric!(grid, 1.0, 1.25)
    column_step!(thermodynamics, insulating_fluxes(), 10.0)

    @test column_values(thermodynamics.fields.internal_energy) ≈ expected_expanded_profile

    basal_grid = RectilinearGrid(size = 4,
                                 z = SeaIceColumnDiscretization((0, 1)),
                                 topology = (Flat, Flat, Bounded))
    # Congelation lays down ice at the basal Dirichlet temperature; the swept enthalpy is the BL99 new-ice
    # enthalpy internal_energy(relation, Tᵇ, S). With zero conductivity the Dirichlet face adds nothing to the
    # interior solve, isolating the conservative basal-fill remap.
    basal_bottom_temperature = -1.8
    fill_energy = internal_energy(relation, basal_bottom_temperature, 0.0)
    basal_growth = prescribed_salinity_enthalpy_thermodynamics(basal_grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = no_transport,
        heat_boundary_conditions = (top = FluxBoundary(),
                                    bottom = PrescribedTemperature(basal_bottom_temperature)))

    set_column_metric!(basal_grid, 1.0, 1.0; surface = 1.0, move = :base)
    set!(basal_growth; bulk_salinity = 0.0, temperature = z -> -12 + 8z)
    basal_initial_profile = column_values(basal_growth.fields.internal_energy)
    basal_growth_faces = collect(range(-0.25, 1.0; length = 5))
    expected_basal_growth_profile =
        conservative_column_remap(basal_initial_profile,
                                  source_faces,
                                  basal_growth_faces;
                                  fill_value = fill_energy)

    set_column_metric!(basal_grid, 1.0, 1.25; surface = 1.0, move = :base)
    column_step!(basal_growth, insulating_fluxes(), 10.0)

    @test column_values(basal_growth.fields.internal_energy) ≈ expected_basal_growth_profile

    flux_forced = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport = no_transport,
        heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary()))

    set_column_metric!(grid, 1.0, 1.0)
    set!(flux_forced; bulk_salinity = 0.0, temperature = -10.0)
    flux_initial_energy = column_integrated_energy(flux_forced)

    set_column_metric!(grid, 1.0, 1.25)
    column_step!(flux_forced, (top = 2.0, bottom = 0.0), 10.0)

    # Upward-positive top flux removes energy: net change is Δt·(−Qᵘ) = −20.
    @test abs(column_integrated_energy(flux_forced) - 1.25 * flux_initial_energy + 20.0) < 1e-8

    salinity_grid = RectilinearGrid(size = 4,
                                    z = SeaIceColumnDiscretization((0, 1)),
                                    topology = (Flat, Flat, Bounded))

    evolving = evolving_salinity_mushy_thermodynamics(salinity_grid;
        relation,
        energy_transport = no_transport,
        salinity_transport = NoSalinityTransport(),
        heat_boundary_conditions)

    set!(evolving; bulk_salinity = 5.0, temperature = -10.0)
    initial_salt = column_integrated_salinity(evolving)
    initial_salinity = column_values(evolving.fields.bulk_salinity)

    set_column_metric!(salinity_grid, 1.0, 1.25)
    column_salinity_time_step!(evolving, 10.0)

    @test column_integrated_salinity(evolving) ≈ 1.25 * initial_salt
    @test column_values(evolving.fields.bulk_salinity) ≈ initial_salinity

    set_column_metric!(salinity_grid, 1.0, 1.0)
    set!(evolving; bulk_salinity = z -> 2 + 6z, temperature = -10.0)
    initial_salinity = column_values(evolving.fields.bulk_salinity)
    expected_salinity =
        conservative_column_remap(initial_salinity, source_faces, contracted_faces)

    set_column_metric!(salinity_grid, 1.0, 0.75)
    column_salinity_time_step!(evolving, 10.0)

    @test column_values(evolving.fields.bulk_salinity) ≈ expected_salinity
end

@testset "Column energy manufactured convergence" begin
    spatial_errors = [spatial_manufactured_solution_error(N) for N in (16, 32, 64)]
    spatial_slopes = log2.(spatial_errors[1:2] ./ spatial_errors[2:3])

    @test all(spatial_errors .< [30, 8, 2])
    @test all(spatial_slopes .> 1.8)

    temporal_errors = [temporal_manufactured_solution_error(dt) for dt in (0.04, 0.02, 0.01)]
    temporal_slopes = log2.(temporal_errors[1:2] ./ temporal_errors[2:3])

    @test all(temporal_errors .< [6000, 3200, 1700])
    @test all(temporal_slopes .> 0.9)
end

salinity_variance(thermodynamics) = begin
    salinity = column_values(thermodynamics.fields.bulk_salinity)
    salinity_mean = sum(salinity) / length(salinity)
    sum(abs2, salinity .- salinity_mean) / length(salinity)
end

@testset "Column salinity fixed-grid time step" begin
    grid = RectilinearGrid(size = 8,
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Flat, Flat, Bounded))

    relation = QuadraticLiquidusEnergyRelation(Float64)
    heat_boundary_conditions = insulating_boundaries()
    energy_transport = ConductiveTemperatureTransport(conductivity = 2.0)
    salinity_transport = BulkSalinityDiffusion(diffusivity = 1e-4)
    dt = 100

    evolving = evolving_salinity_mushy_thermodynamics(grid;
        relation,
        energy_transport,
        salinity_transport,
        heat_boundary_conditions)

    set!(evolving;
         bulk_salinity = z -> 5 + sin(2pi * z),
         temperature = -10.0)

    initial_salt = column_integrated_salinity(evolving)
    initial_salinity = column_values(evolving.fields.bulk_salinity)
    initial_variance = salinity_variance(evolving)

    column_salinity_time_step!(evolving, dt)

    @test abs(column_integrated_salinity(evolving) - initial_salt) < 1e-12
    @test salinity_variance(evolving) <= initial_variance + 1e-13
    @test maximum(abs.(column_values(evolving.fields.bulk_salinity) .- initial_salinity)) > 1e-3

    fixed = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport,
        heat_boundary_conditions)

    no_transport = evolving_salinity_mushy_thermodynamics(grid;
        relation,
        energy_transport,
        salinity_transport = NoSalinityTransport(),
        heat_boundary_conditions)

    set!(fixed;
         bulk_salinity = z -> 4 + z,
         temperature = z -> -15 + 3z)

    set!(no_transport;
         bulk_salinity = z -> 4 + z,
         temperature = z -> -15 + 3z)

    initial_no_transport_salinity = column_values(no_transport.fields.bulk_salinity)

    for _ in 1:100
        column_step!(fixed, insulating_fluxes(), 1000)
        column_step!(no_transport, insulating_fluxes(), 1000)
    end

    @test maximum(abs.(column_values(fixed.fields.internal_energy) .-
                       column_values(no_transport.fields.internal_energy))) < 1e-10
    @test maximum(abs.(column_values(fixed.fields.temperature) .-
                       column_values(no_transport.fields.temperature))) < 1e-10
    @test maximum(abs.(column_values(no_transport.fields.bulk_salinity) .-
                       initial_no_transport_salinity)) == 0
    @test maximum(abs.(column_values(fixed.fields.bulk_salinity) .-
                       initial_no_transport_salinity)) == 0
end

@testset "Steady linear temperature with saline enthalpy" begin
    for FT in (Float32, Float64), salinity in (FT(4), FixedDrainedIceSalinityProfile(FT))
        grid = RectilinearGrid(CPU(), FT; size = 8,
                               z = SeaIceColumnDiscretization((0, 1)),
                               topology = (Flat, Flat, Bounded))
        th = prescribed_salinity_enthalpy_thermodynamics(grid;
            relation = FixedSalinityBrinePocketEnergyRelation(FT),
            salinity_profile = salinity,
            energy_transport = ConductiveTemperatureTransport(conductivity = FT(2)),
            heat_boundary_conditions = (bottom = PrescribedTemperature(FT(-3)),
                                        top = PrescribedTemperature(FT(-0.3))))
        set!(th; bulk_salinity = salinity, temperature = z -> -3 + 2.7z)
        initial_temperature = copy(interior(th.fields.temperature))
        column_step!(th, insulating_fluxes(), FT(100000))
        drift = maximum(abs, interior(th.fields.temperature) .- initial_temperature)
        @test drift <= 100eps(FT) * maximum(abs, initial_temperature)
    end
end

@testset "Column motion across several layers" begin
    for FT in (Float32, Float64), (previous, current) in ((0.1, 0.05), (0.05, 0.1)), faces in ((0, 1), (10, 11))
        grid = RectilinearGrid(CPU(), FT; size = 8, z = SeaIceColumnDiscretization(faces), topology = (Flat, Flat, Bounded))
        relation = FixedSalinityBrinePocketEnergyRelation(FT)
        energy_transport = ConductiveTemperatureTransport(FT; conductivity = 0)
        heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary())
        th = prescribed_salinity_enthalpy_thermodynamics(grid; salinity_profile = 0,
            relation, energy_transport, heat_boundary_conditions)
        initial = FT[-2, -20, -20, -20, -20, -20, -20, -20]
        expected = previous > current ? fill(FT(-20), 8) : FT[-2, -2, -2, -2, -11, -20, -20, -20]
        set!(th; temperature = reshape(initial, 1, 1, 8))
        grid.z.hsⁿ .= 0; grid.z.hs⁻ .= 0
        grid.z.hbⁿ .= FT(-current); grid.z.hb⁻ .= FT(-previous)
        column_energy_time_step!(th, (top = FT(0), bottom = FT(0)), Clock(time = FT(0)), NamedTuple(), ConstantField(FT(0)), FT(5400))
        @test vec(Array(interior(th.fields.temperature))) ≈ expected
        @test Oceananigans.Grids.znode(1, 1, 1, grid, Center(), Center(), Face()) ≈ FT(-current)

        salt = evolving_salinity_mushy_thermodynamics(grid; relation, energy_transport,
            salinity_transport = NoSalinityTransport(), heat_boundary_conditions)
        set!(salt; bulk_salinity = reshape(-initial, 1, 1, 8), temperature = FT(-5))
        column_salinity_time_step!(salt, FT(5400))
        @test vec(Array(interior(salt.fields.bulk_salinity))) ≈ -expected
    end
end
