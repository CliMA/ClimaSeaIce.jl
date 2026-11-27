using Test
using Oceananigans
using ClimaSeaIce

using Oceananigans.Architectures: on_architecture
using Oceananigans.BoundaryConditions: TracerReservoir
using Oceananigans.Grids: Open

@inline exterior_thickness(y, t) = 1 + t

@testset "Open-boundary reservoirs advance with the sea ice clock ($timestepper)" for timestepper in (:ForwardEuler, :SplitRungeKutta3)
    grid = RectilinearGrid(size=(8, 4), x=(0, 1), y=(0, 1), topology=(Open, Bounded, Flat))

    u_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing);
                                    west = NormalFlowBoundaryCondition(0.1), east = NormalFlowBoundaryCondition(0.1))
    u = Field{Face, Center, Nothing}(grid; boundary_conditions = u_bcs)
    v = Field{Center, Face, Nothing}(grid)
    set!(u, 0.1)

    h_bcs = FieldBoundaryConditions(west = ValueBoundaryCondition(exterior_thickness; scheme = TracerReservoir()))
    ℵ_bcs = FieldBoundaryConditions(west = ValueBoundaryCondition(1; scheme = TracerReservoir()))

    model = SeaIceModel(grid; velocities = (; u, v), dynamics = nothing, ice_thermodynamics = nothing,
                        advection = UpwindBiased(order = 1), timestepper,
                        boundary_conditions = (; h = h_bcs, ℵ = ℵ_bcs))

    for _ in 1:4
        time_step!(model, 1)
    end

    # On inflow the reservoir takes the exterior value of the current step
    h = on_architecture(CPU(), model.ice_thickness)
    @test all(h[0, 1:4, 1] .≥ exterior_thickness(0, 3))
    @test all(h[1, 1:4, 1] .> 0)
end
