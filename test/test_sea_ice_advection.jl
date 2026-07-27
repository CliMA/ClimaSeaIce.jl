using ClimaSeaIce
using Oceananigans
using Oceananigans: prognostic_fields
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Operators: Δxᶜᶠᶜ, Δyᶠᶜᶜ

@testset "Sea ice advection" begin
    @info "Running sea ice advection test"

    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))

    model = SeaIceModel(grid, advection=WENO())

    @test architecture(model) == architecture(grid)

    @test !(model.velocities.u isa Nothing)
    @test !(model.velocities.v isa Nothing)

    # test that model runs
    @test begin
        time_step!(model, 1)
        true
    end
end

@testset "Sea ice velocities" begin
    @info "Testing sea ice velocities on a TripolarGrid"
    grid  = TripolarGrid(size = (10, 10, 1))
    model = SeaIceModel(grid)

    @test model.velocities.u.boundary_conditions.north.condition == -1
    @test model.velocities.v.boundary_conditions.north.condition == -1

    grid  = TripolarGrid(size = (10, 10, 1), fold_topology = Oceananigans.Grids.RightFaceFolded)
    model = SeaIceModel(grid)

    @test model.velocities.u.boundary_conditions.north.condition == -1
    @test model.velocities.v.boundary_conditions.north.condition == -1
end

@testset "Sea ice advection with snow" begin
    @info "Running sea ice advection with snow test"

    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))
    snow_thermo = snow_slab_thermodynamics(grid)
    model = SeaIceModel(grid; advection=WENO(), snow_thermodynamics=snow_thermo)

    set!(model, h=1, ℵ=1, hs=0.2)

    @test model.snow_thickness isa Field
    @test :hs in keys(prognostic_fields(model))

    @test begin
        time_step!(model, 1)
        true
    end
end

@testset "Volume-form advection conserves ∫(ℵ·h) dA" begin
    grid = RectilinearGrid(size=(64, 64, 1), x=(0,1), y=(0,1), z=(-1,0), halo=(4,4,4), topology=(Periodic, Periodic, Bounded))

    # Convergent/divergent flow toward (0.5, 0.5)
    u = Field{Face, Center, Center}(grid)
    v = Field{Center, Face, Center}(grid)
    set!(u, (x, y, z) -> -0.1 * sin(2π*x) * cos(2π*y))
    set!(v, (x, y, z) -> -0.1 * cos(2π*x) * sin(2π*y))

    model = SeaIceModel(grid; velocities = (; u, v), advection = WENO(order=7))
    set!(model, h = (x, y) -> 1.0 + 0.5*sin(2π*x)*cos(2π*y),
                ℵ = (x, y) -> 0.5 + 0.3*sin(2π*x)*cos(2π*y))

    V₀ = Field(Integral((model.ice_thickness * model.ice_concentration)))[1, 1, 1]
    for _ in 1:50
        time_step!(model, 0.01)
    end
    V₅₀ = Field(Integral((model.ice_thickness * model.ice_concentration)))[1, 1, 1]

    @test isapprox(V₅₀, V₀; rtol = 1e-12)
end

#####
##### Incremental remapping
#####

remapping_grid(Nx; halo = 4) =
    RectilinearGrid(size = (Nx, Nx, 1), x = (0, 1), y = (0, 1), z = (-1, 0),
                    halo = (halo, halo, halo), topology = (Periodic, Periodic, Bounded))

# A discretely non-divergent swirl: u = -∂ψ/∂y, v = ∂ψ/∂x from a streamfunction at cell corners.
function swirling_velocities(grid; amplitude = 0.1)
    ψ = Field{Face, Face, Center}(grid)
    set!(ψ, (x, y, z) -> amplitude / 2π * sin(2π * x) * sin(2π * y))
    fill_halo_regions!(ψ)

    u = Field{Face, Center, Center}(grid)
    v = Field{Center, Face, Center}(grid)

    Nx, Ny, _ = size(grid)
    for i in 1:Nx, j in 1:Ny
        u[i, j, 1] = - (ψ[i, j+1, 1] - ψ[i, j, 1]) / Δyᶠᶜᶜ(i, j, 1, grid)
        v[i, j, 1] = + (ψ[i+1, j, 1] - ψ[i, j, 1]) / Δxᶜᶠᶜ(i, j, 1, grid)
    end

    fill_halo_regions!(u)
    fill_halo_regions!(v)

    return u, v
end

function remapping_model(grid, u, v; snow = false)
    snow_thermodynamics = snow ? snow_slab_thermodynamics(grid) : nothing
    return SeaIceModel(grid; velocities = (; u, v),
                       advection = IncrementalRemapping(),
                       ice_thermodynamics = nothing,
                       snow_thermodynamics,
                       dynamics = nothing,
                       timestepper = :ForwardEuler)
end

@testset "Incremental remapping requirements" begin
    @info "Testing incremental remapping constructor requirements"

    grid = remapping_grid(8)

    # The transport integrates over the step, so it is a Forward-Euler scheme by construction.
    @test_throws ArgumentError SeaIceModel(grid; advection = IncrementalRemapping())
    @test SeaIceModel(grid; advection = IncrementalRemapping(),
                      timestepper = :ForwardEuler) isa SeaIceModel

    thin = remapping_grid(8; halo = 1)
    @test_throws ArgumentError SeaIceModel(thin; advection = IncrementalRemapping(),
                                           timestepper = :ForwardEuler)
end

@testset "Incremental remapping bounds the recovered thickness" begin
    @info "Testing that incremental remapping keeps h = 𝓋/ℵ bounded"

    grid = remapping_grid(64)
    u, v = swirling_velocities(grid)
    model = remapping_model(grid, u, v; snow = true)

    # Discontinuous concentration carried by a sheared flow, with uniform thickness: the exact answer
    # keeps h ≡ 1.5 and hs ≡ 0.3 wherever ice is present. Independently advecting ℵ and 𝓋 = ℵh cannot
    # do this — the ratio is unbounded as ℵ → 0.
    iced(x, y) = ((x - 1/2)^2 + (y - 1/2)^2) < 0.2^2
    set!(model, ℵ = (x, y) -> iced(x, y) ? 1.0 : 0.0,
                h = (x, y) -> iced(x, y) ? 1.5 : 0.0,
                hs = (x, y) -> iced(x, y) ? 0.3 : 0.0)

    for _ in 1:200
        time_step!(model, 0.02)
    end

    h = Array(interior(model.ice_thickness))
    hs = Array(interior(model.snow_thickness))
    ℵ = Array(interior(model.ice_concentration))

    @test maximum(h) ≤ 1.5 * (1 + 1e-4)
    @test maximum(hs) ≤ 0.3 * (1 + 1e-4)
    @test maximum(ℵ) ≤ 1
    @test minimum(ℵ) ≥ 0
end

@testset "Incremental remapping is exact for uniform translation" begin
    @info "Testing incremental remapping under uniform translation"

    grid = remapping_grid(48)

    u = Field{Face, Center, Center}(grid); set!(u, 0.1); fill_halo_regions!(u)
    v = Field{Center, Face, Center}(grid); set!(v, 0.037); fill_halo_regions!(v)

    model = remapping_model(grid, u, v)

    iced(x, y) = ((x - 1/2)^2 + (y - 1/2)^2) < 0.2^2
    set!(model, ℵ = (x, y) -> iced(x, y) ? 1.0 : 0.0,
                h = (x, y) -> iced(x, y) ? 1.5 : 0.0)

    for _ in 1:100
        time_step!(model, 0.05)
    end

    # Without shear the swept region never straddles a cell boundary and the transport is exact, so
    # the thickness is bounded to round-off.
    @test maximum(interior(model.ice_thickness)) ≤ 1.5 * (1 + 1e-12)
end

@testset "Incremental remapping conserves ∫(ℵ·h) dA" begin
    @info "Testing incremental remapping conservation"

    grid = remapping_grid(48)
    u, v = swirling_velocities(grid)
    model = remapping_model(grid, u, v; snow = true)

    set!(model, ℵ = (x, y) -> 0.5 + 0.4 * sin(2π * x) * cos(2π * y), h = 1.5, hs = 0.3)

    content(m) = Field(Integral(m.ice_concentration * m.ice_thickness))[1, 1, 1]
    snow(m) = Field(Integral(m.ice_concentration * m.snow_thickness))[1, 1, 1]

    𝓋₀, 𝓋s₀ = content(model), snow(model)

    for _ in 1:100
        time_step!(model, 0.02)
    end

    @test isapprox(content(model), 𝓋₀; rtol = 1e-12)
    @test isapprox(snow(model), 𝓋s₀; rtol = 1e-12)

    # A uniform thickness stays uniform: the transported content and area come from the same region.
    @test maximum(interior(model.ice_thickness)) ≈ 1.5
    @test minimum(interior(model.ice_thickness)) ≈ 1.5
end

@testset "Incremental remapping needs no minimum concentration" begin
    @info "Testing that incremental remapping conserves content without emptying cells"

    # Flux-form advection empties any cell whose concentration collapses, because the thickness it
    # recovers there is meaningless. That discards content. Incremental remapping bounds the ratio at
    # any concentration, so it keeps the floor at zero and conserves to round-off even when the
    # concentration is discontinuous.
    @test ClimaSeaIce.minimum_ice_concentration(Float64, IncrementalRemapping()) == 0
    @test ClimaSeaIce.minimum_ice_concentration(Float64, WENO(order=5)) > 0

    grid = remapping_grid(48)
    u, v = swirling_velocities(grid)
    model = remapping_model(grid, u, v)

    iced(x, y) = ((x - 1/2)^2 + (y - 1/2)^2) < 0.2^2
    set!(model, ℵ = (x, y) -> iced(x, y) ? 1.0 : 0.0,
                h = (x, y) -> iced(x, y) ? 1.5 : 0.0)

    content(m) = Field(Integral(m.ice_concentration * m.ice_thickness))[1, 1, 1]
    𝓋₀ = content(model)

    for _ in 1:200
        time_step!(model, 0.02)
    end

    @test isapprox(content(model), 𝓋₀; rtol = 1e-12)
    @test minimum(interior(model.ice_concentration)) ≥ 0
end

@testset "Incremental remapping is second-order accurate" begin
    @info "Testing incremental remapping convergence rate"

    function translation_error(Nx; courant = 0.25)
        grid = remapping_grid(Nx)

        u = Field{Face, Center, Center}(grid); set!(u, 1); fill_halo_regions!(u)
        v = Field{Center, Face, Center}(grid); set!(v, 0); fill_halo_regions!(v)

        model = remapping_model(grid, u, v)
        set!(model, ℵ = (x, y) -> 0.5 + 0.4 * sin(2π * x) * cos(2π * y), h = 1.5)

        ℵ₀ = Array(interior(model.ice_concentration))

        steps = round(Int, Nx / courant)
        for _ in 1:steps
            time_step!(model, 1 / steps)
        end

        ℵ = Array(interior(model.ice_concentration))
        return sum(abs, ℵ .- ℵ₀) / length(ℵ)
    end

    errors = translation_error.((16, 32, 64))
    orders = log2.(errors[1:end-1] ./ errors[2:end])

    @test all(orders .> 1.8)
end

@testset "Incremental remapping across an immersed boundary" begin
    @info "Testing incremental remapping against an immersed island"

    underlying_grid = remapping_grid(48)
    island(x, y, z) = ((x - 1/2)^2 + (y - 1/2)^2) < 0.12^2
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBoundary(island))

    u = Field{Face, Center, Center}(grid); set!(u, 0.1); fill_halo_regions!(u)
    v = Field{Center, Face, Center}(grid); set!(v, 0.037); fill_halo_regions!(v)

    model = remapping_model(grid, u, v)
    set!(model, ℵ = 0.3, h = 1.5)
    Oceananigans.TimeSteppers.update_state!(model)

    content(m) = sum(interior(m.ice_concentration, :, :, 1) .* interior(m.ice_thickness, :, :, 1))
    𝓋₀ = content(model)

    for _ in 1:200
        time_step!(model, 0.02)
    end

    ℵ = Array(interior(model.ice_concentration, :, :, 1))
    solid = [ClimaSeaIce.submerged(i, j, 1, grid) for i in 1:size(grid, 1), j in 1:size(grid, 2)]

    # Nothing crosses an immersed face and no reconstruction reaches into land, so ice driven onto a
    # coast is conserved rather than deposited on it.
    @test maximum(ℵ[solid]; init = 0.0) == 0
    @test isapprox(content(model), 𝓋₀; rtol = 1e-10)
    @test minimum(ℵ) ≥ 0
    @test all(isfinite, interior(model.ice_thickness))
end

@testset "Incremental remapping on a TripolarGrid" begin
    @info "Testing incremental remapping across the tripolar fold"

    grid = TripolarGrid(size = (40, 40, 1), z = (-1, 0), halo = (4, 4, 4))

    bcs = name -> ClimaSeaIce.default_sea_ice_boundary_conditions(grid, name)
    u = Field{Face, Center, Nothing}(grid, boundary_conditions = bcs(:u))
    v = Field{Center, Face, Nothing}(grid, boundary_conditions = bcs(:v))
    set!(u, (λ, φ) -> 0.2 * sind(2λ) * cosd(φ))
    set!(v, (λ, φ) -> 0.1 * cosd(2λ))
    fill_halo_regions!(u)
    fill_halo_regions!(v)

    model = SeaIceModel(grid; velocities = (; u, v),
                        advection = IncrementalRemapping(),
                        ice_thermodynamics = nothing,
                        dynamics = nothing,
                        timestepper = :ForwardEuler)

    iced(φ) = φ > 55
    set!(model, ℵ = (λ, φ) -> iced(φ) ? 0.3 : 0.0, h = (λ, φ) -> iced(φ) ? 1.5 : 0.0)

    for _ in 1:100
        time_step!(model, 3600)
    end

    h = Array(interior(model.ice_thickness))
    ℵ = Array(interior(model.ice_concentration))

    # The reconstruction gradients are computed directly in the halo, which gives them the sign the
    # fold requires without a separate signed halo exchange.
    @test maximum(h) ≤ 1.5 * (1 + 1e-8)
    @test maximum(ℵ) ≤ 1
    @test minimum(ℵ) ≥ 0
    @test all(isfinite, h)
end

@testset "Sea ice momentum equations" begin
    @info "Running sea ice momentum equations test"

    grid = RectilinearGrid(size=(10, 10), x=(0, 1), y=(0, 1), topology=(Bounded, Bounded, Flat))
    dynamics = SeaIceMomentumEquation(grid, rheology=ViscousRheology(ν=1000))

    model = SeaIceModel(grid; dynamics, ice_thermodynamics=nothing, advection=WENO())

    @test !(model.velocities.u isa Nothing)
    @test !(model.velocities.v isa Nothing)

    # test that model runs
    @test begin
        time_step!(model, 1)
        true
    end
end
