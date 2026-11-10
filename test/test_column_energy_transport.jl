using ClimaSeaIce.SeaIceThermodynamics:
    ColumnEnergyThermodynamics,
    SeaIceColumnDiscretization,
    ConductiveTemperatureTransport,
    FluxBoundary,
    NoSalinityTransport,
    PrescribedBulkSalinity,
    PrognosticBulkSalinity,
    evolving_salinity_mushy_thermodynamics,
    internal_energy,
    temperature
using ClimaSeaIce: ClimaSeaIce, SeaIceModel
using ClimaSeaIce: IncrementalRemapping
using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.TimeSteppers: time_step!, update_state!
using Oceananigans.Fields: interior, set!
using Oceananigans.Operators: Azᶜᶜᶜ, Δyᶜᶠᶜ
using Test

@testset "Insulated columns conserve heat across the tripolar fold" begin
    for FT in (Float32, Float64), salinity_closure in (PrescribedBulkSalinity(4), PrognosticBulkSalinity())
        grid = TripolarGrid(CPU(), FT; size = (8, 8, 2), halo = (2, 2, 2),
            z = SeaIceColumnDiscretization((0, 1)), fold_topology = Oceananigans.Grids.RightFaceFolded)
        ice = ColumnEnergyThermodynamics(grid; salinity_closure,
            energy_transport = ConductiveTemperatureTransport(conductivity = 2),
            heat_boundary_conditions = insulating_boundaries())
        model = SeaIceModel(grid; ice_thermodynamics = ice, timestepper = :ForwardEuler)
        set!(model; h = 1, ℵ = 1)
        set!(ice; bulk_salinity = 4, temperature = (longitude, latitude, z) -> -5 + z)
        E = ice.fields.internal_energy
        initial_energy = Array(interior(E))
        Nx, Ny, Nz = size(grid)
        Hx, Hy, Hz = Oceananigans.Grids.halo_size(grid)
        fold = (Hx+1:Hx+Nx, Hy+Ny+1, Hz+1:Hz+Nz)
        mirror = (Hx+Nx:-1:Hx+1, Hy+Ny, Hz+1:Hz+Nz)
        time_step!(model, 60)
        temperature = Array(parent(ice.fields.temperature))
        @test sum(Array(interior(E)); dims = 3) ≈ sum(initial_energy; dims = 3)
        @test temperature[fold...] ≈ temperature[mirror...]
    end
end

@testset "Column properties follow conservative ice-volume remapping" begin
    N = 16
    grid = RectilinearGrid(Float64; size = (N, N, 3), x = (0, 1), y = (0, 1),
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Periodic, Periodic, Bounded), halo = (4, 4, 4))
    thermodynamics = evolving_salinity_mushy_thermodynamics(grid;
        energy_transport = ConductiveTemperatureTransport(conductivity = 0.0),
        salinity_transport = NoSalinityTransport(),
        heat_boundary_conditions = (top = FluxBoundary(), bottom = FluxBoundary()))
    model = SeaIceModel(grid; ice_thermodynamics = thermodynamics,
                        advection = IncrementalRemapping(), timestepper = :ForwardEuler,
                        top_heat_flux = 0.0, bottom_heat_flux = 0.0)
    anomaly(x, y) = exp(-((x - 0.35)^2 + (y - 0.45)^2) / 0.1^2)
    set!(model; h = 1.0, ℵ = (x, y) -> 0.5anomaly(x, y),
         u = (x, y) -> 0.2, v = (x, y) -> 0.1)
    set!(thermodynamics; bulk_salinity = (x, y, z) -> 4 + anomaly(x, y),
         temperature = (x, y, z) -> -5 - anomaly(x, y))

    E⁰ = Array(interior(thermodynamics.fields.internal_energy))
    V⁰ = Array(interior(model.ice_thickness)) .* Array(interior(model.ice_concentration))
    S⁰ = Array(interior(thermodynamics.fields.bulk_salinity))
    energy⁰ = sum(E⁰ .* V⁰)
    salt⁰ = sum(S⁰ .* V⁰)

    for _ in 1:20
        time_step!(model, 0.05)
    end

    E = Array(interior(thermodynamics.fields.internal_energy))
    S = Array(interior(thermodynamics.fields.bulk_salinity))
    V = Array(interior(model.ice_thickness)) .* Array(interior(model.ice_concentration))
    @test isapprox(sum(V), sum(V⁰); rtol = 1e-12)
    @test isapprox(sum(E .* V), energy⁰; rtol = 1e-12)
    @test isapprox(sum(S .* V), salt⁰; rtol = 1e-12)
    @test maximum(abs.(E .- E⁰)) > 1e-6
    @test all(isfinite, E) && all(isfinite, S)
end

@testset "Two-column transport carries energy and salt with ice volume" begin
    grid = RectilinearGrid(size = (2, 2, 2), x = (0, 2), y = (0, 2),
                           z = SeaIceColumnDiscretization((0, 1)),
                           topology = (Bounded, Bounded, Bounded), halo = (2, 2, 2))
    ice = evolving_salinity_mushy_thermodynamics(grid)
    model = SeaIceModel(grid; ice_thermodynamics = ice, advection = IncrementalRemapping(), timestepper = :ForwardEuler)
    set!(model; h = 1, ℵ = 1)
    set!(ice; internal_energy = (x, y, z) -> ifelse(x < 1, -1e7, -2e7),
              bulk_salinity = (x, y, z) -> ifelse(x < 1, 4, 0))
    transports = model.advection.transports
    set!(transports.x.𝓋, (x, y) -> ifelse(x == 1, 0.25, 0))
    set!(model.timestepper.Gⁿ.h, (x, y) -> ifelse(x < 1, -0.25, 0.25))
    set!(model.timestepper.Gⁿ.ℵ, (x, y) -> ifelse(x < 1, -0.25, 0.25))
    ClimaSeaIce.compute_column_transports!(ice, CPU(), grid, transports)
    ClimaSeaIce.dynamic_time_step!(model, 1)
    V = Array(interior(model.ice_thickness)) .* Array(interior(model.ice_concentration))
    E = Array(interior(ice.fields.internal_energy))
    S = Array(interior(ice.fields.bulk_salinity))
    @test vec(V) ≈ [0.75, 1.25, 0.75, 1.25]
    @test E[:, 1, 1] ≈ [-1e7, -1.8e7]
    @test S[:, 1, 1] ≈ [4, 0.8]
    @test isapprox(sum(E .* V) / 4, -3e7; rtol = 1e-14)
    @test isapprox(sum(S .* V) / 4, 4; rtol = 1e-14)
end

@testset "Layered remapping through a tripolar fold conserves volume and heat" begin
    # Thin ice facing thick ice across the fold: the state that took one column of a global run below absolute zero.
    function remap_facing_columns(direction; steps = 6, Δt = 5400.0)
        grid = TripolarGrid(CPU(), Float64; size = (8, 8, 8), halo = (3, 3, 3),
                            z = SeaIceColumnDiscretization((0, 1)), fold_topology = Oceananigans.Grids.RightFaceFolded)
        ice = ColumnEnergyThermodynamics(grid; salinity_closure = PrescribedBulkSalinity(4),
                                         energy_transport = ConductiveTemperatureTransport(conductivity = 2),
                                         heat_boundary_conditions = insulating_boundaries())
        model = SeaIceModel(grid; ice_thermodynamics = ice, timestepper = :ForwardEuler,
                            advection = IncrementalRemapping(), dynamics = nothing)
        Nx, Ny, Nz = size(grid)
        left, right = 3, Nx - 2
        donor, receiver = direction === :to_right ? (left, right) : (right, left)
        h, ℵ, v = model.ice_thickness, model.ice_concentration, model.velocities.v
        E, T = ice.fields.internal_energy, ice.fields.temperature

        interior(h) .= 0
        interior(ℵ) .= 0
        h[donor, Ny, 1], ℵ[donor, Ny, 1] = 1.46, 0.626
        h[receiver, Ny, 1], ℵ[receiver, Ny, 1] = 1.37, 0.0118
        set!(ice; bulk_salinity = 4, temperature = (λ, φ, z) -> -1)
        for k in 1:Nz
            T[receiver, Ny, k] = -3
        end
        ClimaSeaIce.SeaIceThermodynamics.compute_column_internal_energy!(ice)

        # Moves about 0.7% of the donor's area per step through the fold; the left half of the fold row is the one that is kept.
        speed = 0.00706 * Δyᶜᶠᶜ(left, Ny + 1, 1, grid) / Δt
        v[left, Ny + 1, 1] = direction === :to_right ? speed : -speed
        fill_halo_regions!(v)
        update_state!(model)

        volume(i) = h[i, Ny, 1] * ℵ[i, Ny, 1]
        total_volume() = sum(volume(i) * Azᶜᶜᶜ(i, Ny, 1, grid) for i in 1:Nx)
        total_heat() = sum(volume(i) * Azᶜᶜᶜ(i, Ny, 1, grid) * sum(E[i, Ny, k] for k in 1:Nz) / Nz for i in 1:Nx)
        layer_energy() = [E[i, Ny, k] for i in (donor, receiver), k in 1:Nz if volume(i) > 0]

        V₀, Q₀ = total_volume(), total_heat()
        E₋, E₊ = extrema(layer_energy())
        receiver_volume₀, donor_volume₀ = volume(receiver), volume(donor)
        energy_in_range = true
        for _ in 1:steps
            time_step!(model, Δt)
            energy_in_range &= all(e -> E₋ * (1 - 1e-12) ≤ e ≤ E₊ * (1 + 1e-12), layer_energy())
        end

        return (; transports = model.advection.transports, V₀, Q₀, ΔV = total_volume() - V₀, ΔQ = total_heat() - Q₀, energy_in_range,
                gained = volume(receiver) - receiver_volume₀, lost = donor_volume₀ - volume(donor),
                receiver_volume = volume(receiver), receiver_energy = E[receiver, Ny, 1])
    end

    to_right = remap_facing_columns(:to_right)
    to_left = remap_facing_columns(:to_left)

    for direction in (to_right, to_left)
        @test abs(direction.ΔV) ≤ 1e-12 * direction.V₀
        @test abs(direction.ΔQ) ≤ 1e-12 * abs(direction.Q₀)
        @test direction.energy_in_range
        @test direction.gained > 0
        @test direction.lost > 0
    end

    # A mirror image of the flow gives a mirror image of the result.
    @test to_right.receiver_volume ≈ to_left.receiver_volume rtol = 1e-12
    @test to_right.receiver_energy ≈ to_left.receiver_energy rtol = 1e-12

    for transports in to_right.transports, name in (:ℵ, :𝓋, :𝓋s, :E, :bulk_salinity)
        @test transports[name].boundary_conditions.north.condition == -1
    end
end
