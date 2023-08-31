#####
##### Surface thermal boundary conditions
#####

struct MeltingConstrainedFluxBalance{STS}
    surface_temperature_solver :: STS
end

struct LinearizedSurfaceTemperatureSolver end
struct NonlinearSurfaceTemperatureSolver end

"""
    MeltingConstrainedFluxBalance(surface_temperature_solver=NonlinearSurfaceTemperatureSolver())

Return a boundary condition that determines the surface temperature ``T̃``
to equilibrate the surface heat fluxes,

```math
Q₁(Tᵤ) + Q₂(Tᵤ) + ⋯ - Qᵢ = Σᴺ Qₙ(Tᵤ) - Qᵢ = δQ ,
```

where ``Qᵢ`` is the intrinsic flux into the surface from within the ice
(typically, a conductive flux), ``Qₙ`` represent external fluxes into the
air above the ice, ``δQ`` is the residual flux, and ``Tᵤ`` is the (upper)
surface temperature. ``Tᵤ`` is evaluated under the constraint that

```math
Tᵤ ≤ Tₘ(S),
```

where ``Tₘ(S)`` is the melting temperature, which is a function of the
ice salinity at the surface, ``S``. When ``Tᵤ < Tₘ(S)``, the surface is
frozen and ``δQ = 0``. When the constraint operates, such that ``Tᵤ = Tₘ(S)``,
the surface is melting and the residual flux is non-zero.

```math
δQ ≡ Σᴺ Qₙ(Tₘ) - Qᵢ(Tₘ).
```

The residual flux is consumed by the cost of transforming ice into liquid water,
and is related to the rate of change of ice thickness, ``h``, by

```math
d/dt hᵤ = δQ / ℒ(Tᵤ)
```

where ``ℒ(Tᵤ)`` is the latent heat, equal to the different between
the higher internal energy of liquid water and the lower internal energy of solid ice,
at the temperature ``Tᵤ``.

"""
MeltingConstrainedFluxBalance() = MeltingConstrainedFluxBalance(NonlinearSurfaceTemperatureSolver())

#####
##### Flux imbalance and temperature
#####

@inline function top_flux_imbalance(i, j, grid, top_thermal_bc, top_temperature,
                                    internal_fluxes, external_fluxes, clock, model_fields)

    # Schematic of the Stefan condition
    #        
    #  air       ↑   Qx ≡ external_fluxes. Example: Qx = σ T⁴
    #          |⎴⎴⎴|
    # ----------------------- ↕ hᵤ(t) → ℒ ∂t hᵤ = δQ | T ≤ Tₘ
    #          |⎵⎵⎵|
    #  ice       ↑   Qi ≡ internal_fluxes. Example Qi = - k ∂z T
    #      

    Qi = getflux(internal_fluxes, i, j, grid, top_temperature, clock, model_fields)
    Qx = getflux(external_fluxes, i, j, grid, top_temperature, clock, model_fields)

    # The imbalance is defined such that
    # negative imbalance => heat accumulates (out > in) ⟹  growth
    return Qx - Qi
end

@inline top_temperature(i, j, grid, ::PrescribedTemperature, Tu, args...) = Tu

using RootSolvers: SecantMethod, find_zero, CompactSolution

@inline function top_temperature(i, j, grid, top_thermal_bc,
                                 current_surface_temperature,
                                 internal_fluxes, external_fluxes,
                                 clock, model_fields)

    Tu = @inbounds current_surface_temperature[i, j, 1]
    T₁ = Tu
    T₂ = Tu - 1
    FT = eltype(grid)
    method = SecantMethod{FT}(T₁, T₂)
    solution_type = CompactSolution()

    flux_balance(T) = getflux(external_fluxes, i, j, grid, T, clock, model_fields) -
                      getflux(internal_fluxes, i, j, grid, T, clock, model_fields)
                      
    solution = find_zero(flux_balance, method, solution_type)
    
    return solution.root
end



