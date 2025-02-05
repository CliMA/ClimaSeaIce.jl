using Oceananigans
using ClimaSeaIce
using ClimaSeaIce.Rheologies: CenteredWENO5, interpolate_xᶠ, interpolate_yᶠ
using Oceananigans.Advection: AbstractCenteredAdvectionScheme, LeftBias, RightBias

using Oceananigans
using Oceananigans.Models.ShallowWaterModels: VectorInvariantFormulation, ConservativeFormulation
using JLD2

import Oceananigans.Advection: _symmetric_interpolate_xᶠᵃᵃ, 
                               _symmetric_interpolate_yᵃᶠᵃ, 
                               _symmetric_interpolate_zᵃᵃᶠ,
                               _biased_interpolate_xᶠᵃᵃ

struct MyCenteredWENO5  <: AbstractCenteredAdvectionScheme{3, Float64} end
struct MyCenteredWENO52 <: AbstractCenteredAdvectionScheme{3, Float64} end

@inline _symmetric_interpolate_xᶠᵃᵃ(i, j, k, grid, ::MyCenteredWENO5, c) = interpolate_xᶠ(i, j, k, grid, CenteredWENO5(), c)
@inline _symmetric_interpolate_yᵃᶠᵃ(i, j, k, grid, ::MyCenteredWENO5, c) = zero(grid)
@inline _symmetric_interpolate_zᵃᵃᶠ(i, j, k, grid, ::MyCenteredWENO5, c) = zero(grid)

@inline _symmetric_interpolate_xᶠᵃᵃ(i, j, k, grid, ::MyCenteredWENO52, c) = 
    (_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, WENO(), LeftBias(), c) + _biased_interpolate_xᶠᵃᵃ(i, j, k, grid, WENO(), RightBias(), c)) / 2
@inline _symmetric_interpolate_yᵃᶠᵃ(i, j, k, grid, ::MyCenteredWENO52, c) = zero(grid)
@inline _symmetric_interpolate_zᵃᵃᶠ(i, j, k, grid, ::MyCenteredWENO52, c) = zero(grid)


"""
This simulation is a simple 1D advection to test the 
validity of the stretched WENO scheme
"""

arch = CPU()

#parameters
Nx = 200

# 1D grid constructions
grid = RectilinearGrid(size = Nx, x = (-1, 1), halo = 7, topology = (Periodic, Flat, Flat))

# the initial condition
@inline G(x, β, z) = exp(-β*(x - z)^2)
@inline F(x, α, a) = √(max(1 - α^2*(x-a)^2, 0.0))

Z = -0.7
δ = 0.005
β = log(2)/(36*δ^2)
a = 0.5
α = 10

@inline function c₀_1D(x) 
    if x <= -0.6 && x >= -0.8
        return 1/6*(G(x, β, Z-δ) + 4*G(x, β, Z) + G(x, β, Z+δ))
    elseif x <= -0.2 && x >= -0.4
        return 1.0
    elseif x <= 0.2 && x >= 0.0
        return 1.0 - abs(10 * (x - 0.1))
    elseif x <= 0.6 && x >= 0.4
        return 1/6*(F(x, α, a-δ) + 4*F(x, α, a) + F(x, α, a+δ))
    else
        return 0.0
    end
end

# c₀_1D(x) = sin(4π * xs

Δt_max   = 0.2 * minimum_xspacing(grid)
end_time = 2.0

tot_iter = end_time ÷ Δt_max

c_real = CenterField(grid)
    
model = ShallowWaterModel(; grid = grid,
                         tracers = :c,
              momentum_advection = nothing,
                  mass_advection = nothing,
                tracer_advection = MyCenteredWENO5(), #, # WENO(order=5),
      gravitational_acceleration = 1.0,
                     formulation = VectorInvariantFormulation())

set!(model, h=1.0, u=-1.0)
set!(model, c=c₀_1D)

c  = model.tracers.c
cs = []

for i = 1:tot_iter+1
    time_step!(model, Δt_max)
    push!(cs, deepcopy(interior(c, :, 1, 1)))
end
