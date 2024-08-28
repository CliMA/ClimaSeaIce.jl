using Oceananigans
using Oceananigans.Units
using ClimaSeaIce

top_temperature = -10 # C
conductivity = 2 # W m² K⁻¹ .. ?

grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

top_heat_boundary_condition = PrescribedTemperature(top_temperature)
internal_heat_flux = ConductiveFlux(; conductivity)

model = SlabSeaIceModel(grid; top_heat_boundary_condition, internal_heat_flux)

