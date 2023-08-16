using Test

using ClimaSeaIce
using ClimaSeaIce.EulerianThermodynamicSeaIceModels: EulerianThermodynamicSeaIceModel, MolecularDiffusivity

κ = 1e-5
grid = RectilinearGrid(size=3, z=(-3, 0), topology=(Flat, Flat, Bounded))
closure = MolecularDiffusivity(grid, κ_ice=κ, κ_water=κ)
model = EulerianThermodynamicSeaIceModel(; grid, closure)

@test model isa EulerianThermodynamicSeaIceModel

# Test that it runs
simulation = Simulation(model; Δt = 0.1 / κ)

@test begin
    try
        run!(simulation)
        true
    catch
        false
    end
end

