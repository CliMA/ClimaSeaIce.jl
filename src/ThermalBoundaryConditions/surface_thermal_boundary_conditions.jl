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

Return a boundary condition that determines the surface temperature ``Tₛ``
to equilibrate the surface heat fluxes,

```math
Qₓ₁(Tₛ) + Qₓ₂(Tₛ) + ⋯ - Qᵢ = Σᴺ Qₓₙ(Tₛ) - Qᵢ = δQ ,
```

where ``Qᵢ`` is the intrinsic flux into the surface from within the ice
(typically, a conductive flux), ``Qₓₙ`` represent external fluxes into the
air above the ice, ``δQ`` is the residual flux, and ``Tₛ`` is the (upper)
surface temperature. ``Tₛ`` is evaluated under the constraint that

```math
Tₛ ≤ Tₘ(S),
```

where ``Tₘ(S)`` is the melting temperature, which is a function of the
ice salinity at the surface, ``S``. When ``Tₛ < Tₘ(S)``, the surface is
frozen and ``δQ = 0``. When the constraint operates, such that ``Tₛ = Tₘ(S)``,
the surface is melting and the residual flux is non-zero.

```math
δQ ≡ Σᴺ Qₙ(Tₘ) - Qᵢ(Tₘ).
```

The residual flux is consumed by the cost of transforming ice into liquid water,
and is related to the rate of change of ice thickness, ``h``, by

```math
d/dt hₛ = δQ / ℒ(Tₛ)
```

where ``ℒ(Tₛ)`` is the latent heat, equal to the different between
the higher internal energy of liquid water and the lower internal energy of solid ice,
at the temperature ``Tₛ``.

"""
MeltingConstrainedFluxBalance() = MeltingConstrainedFluxBalance(NonlinearSurfaceTemperatureSolver())

#####
##### Flux imbalance and temperature
#####

@inline function surface_flux_imbalance(i, j, grid, surface_thermal_bc, surface_temperature,
                                        internal_fluxes, external_fluxes, clock, model_fields)

    # Schematic of the Stefan condition at an upper surface
    #        
    #  air       ↑   Qx ≡ external_fluxes. Example: Qx = σ T⁴
    #          |⎴⎴⎴|
    # ----------------------- ↕ hₛ(t) → ℒ ∂t hₛ = δQ | T ≤ Tₘ
    #          |⎵⎵⎵|
    #  ice       ↑   Qi ≡ internal_fluxes. Example Qi = - k ∂z T
    #      

    Qi = getflux(internal_fluxes, i, j, grid, surface_temperature, clock, model_fields)
    Qx = getflux(external_fluxes, i, j, grid, surface_temperature, clock, model_fields)

    # The imbalance is defined such that
    # negative imbalance => heat accumulates (out > in) ⟹  growth
    return Qx - Qi
end

@inline surface_temperature(i, j, grid, ::PrescribedTemperature, Tu, args...) = Tu

using RootSolvers: SecantMethod, find_zero, CompactSolution

@inline function surface_temperature(i, j, grid, surface_thermal_bc,
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



