using ClimaSeaIce
using Oceananigans
using Oceananigans: prognostic_fields
using Oceananigans.Architectures: architecture

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
