using Test
using JLD2
using MPI
using Oceananigans
using ClimaSeaIce

const script_path = abspath(joinpath(@__DIR__, "_halo_extension_mpi_worker.jl"))
const result_path = abspath(joinpath(@__DIR__, "_halo_extension_result.jld2"))

test_script = """
using MPI
MPI.Init()

using JLD2
using Test
using Oceananigans
using Oceananigans.DistributedComputations: Distributed, Partition, @root
using Oceananigans.Grids: halo_size, topology
using Oceananigans.Utils: KernelParameters
using ClimaSeaIce
using ClimaSeaIce.SeaIceDynamics: SeaIceMomentumEquation, SplitExplicitSolver
using ClimaSeaIce.Rheologies: ElastoViscoPlasticRheology

arch = Distributed(CPU())
grid = RectilinearGrid(arch; size=(32, 8, 4), extent=(1, 1, 1),
                       halo=(4, 4, 4),
                       topology=(Periodic, Bounded, Bounded))

substeps = 8
solver   = SplitExplicitSolver(grid; substeps)
dynamics = SeaIceMomentumEquation(grid; solver, rheology=ElastoViscoPlasticRheology(eltype(grid)))
model    = ClimaSeaIce.SeaIceModel(grid; dynamics)

Hx, Hy, Hz = halo_size(grid)
Hx = max(substeps + 2, Hx)
@test halo_size(model.velocities.u.grid) == (Hx, Hy, Hz)
@test halo_size(model.dynamics.auxiliaries.fields.P.grid) == (Hx, Hy, Hz)
# Prognostic tracer-like fields must share the extended halos with the EVP
# auxiliaries — otherwise the stress kernel reads h/ℵ out-of-bounds.
@test halo_size(model.ice_thickness.grid)     == (Hx, Hy, Hz)
@test halo_size(model.ice_concentration.grid) == (Hx, Hy, Hz)

Nx, Ny, Nz = size(grid)

kp = model.dynamics.solver.kernel_parameters
# `split_explicit_kernel_size` for a ConnectedTopology in x is `-H+2:N+H-1`;
# in y this run is Bounded, so the kernel just iterates `1:N`.
@test kp == KernelParameters(-Hx+2:Nx+Hx-1, 1:Ny)
"""

@testset "Distributed velocity-halo extension (regression test)" begin
    write(script_path, test_script)
    run(`$(mpiexec()) -n 2 $(Base.julia_cmd()) --project -O0 $(script_path)`)
    rm(script_path)
end
