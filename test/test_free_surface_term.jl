using Test
using Oceananigans
using ClimaSeaIce

using Oceananigans.Fields: interior

# ∂u/∂t = - g ∇η: a surface tilted by a constant slope accelerates ice of any thickness down-slope at
# - g ∇η, so after a forward-Euler step from rest the velocity is exactly - g ∇η Δt.
@testset "Ocean surface tilt accelerates the ice down-slope" begin
    @info "Testing the ocean surface tilt for explicit and split-explicit solvers..."

    grid = RectilinearGrid(size = (8, 8), x = (0, 8e4), y = (0, 8e4), topology = (Bounded, Bounded, Flat))

    ∂η∂x, ∂η∂y = 1e-5, -2e-5
    η = Field{Center, Center, Nothing}(grid)
    set!(η, (x, y) -> ∂η∂x * x + ∂η∂y * y)

    g = 9.5   # not the default, so the keyword is exercised too
    Δt = 60

    for solver in (ExplicitSolver(), SplitExplicitSolver(grid; substeps = 10))
        dynamics = SeaIceMomentumEquation(grid; coriolis = nothing, rheology = ViscousRheology(ν = 0), solver,
                                          ocean_surface_height = η, gravitational_acceleration = g)

        model = SeaIceModel(grid; dynamics, advection = nothing, ice_thermodynamics = nothing,
                            timestepper = :ForwardEuler)
        set!(model, h = 1, ℵ = 1)
        time_step!(model, Δt)

        u = Array(interior(model.velocities.u))
        v = Array(interior(model.velocities.v))

        @test all(u[2:8, :, 1] .≈ - g * ∂η∂x * Δt)
        @test all(v[:, 2:8, 1] .≈ - g * ∂η∂y * Δt)
    end
end
