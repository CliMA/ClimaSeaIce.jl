#####
##### Surface heat boundary conditions
#####

struct MeltingConstrainedFluxBalance{STS}
    top_surface_temperature_solver :: STS
end

struct LinearizedSurfaceTemperatureSolver end
struct NonlinearSurfaceTemperatureSolver end

"""
    MeltingConstrainedFluxBalance(top_surface_temperature_solver=NonlinearSurfaceTemperatureSolver())

Return a boundary condition that determines the top (or "upper surface") temperature ``Tᵤ``
to equilibrate the top heat fluxes,

```math
Qₓ₁(Tᵤ) + Qₓ₂(Tᵤ) + ⋯ - Qᵢ = Σᴺ Qₓₙ(Tᵤ) - Qᵢ = δQ ,
```

where ``Qᵢ`` is the intrinsic flux into the top from within the ice
(typically, a conductive flux), ``Qₓₙ`` represent external fluxes into the
air above the ice, ``δQ`` is the residual flux, and ``Tᵤ`` is the (upper)
top temperature. ``Tᵤ`` is evaluated under the constraint that

```math
Tᵤ ≤ Tₘ(S),
```

where ``Tₘ(S)`` is the melting temperature, which is a function of the
ice salinity at the top, ``S``. When ``Tᵤ < Tₘ(S)``, the top is
frozen and ``δQ = 0``. When the constraint operates, such that ``Tᵤ = Tₘ(S)``,
the top is melting and the residual flux is non-zero.

```math
δQ ≡ Σᴺ Qₙ(Tₘ) - Qᵢ(Tₘ).
```

The residual flux is consumed by the cost of transforming ice into liquid water,
and is related to the rate of change of ice thickness, ``h``, by

```math
\\frac{dhₛ}{dt} = δQ / ℒ(Tᵤ)
```

where ``ℒ(Tᵤ)`` is the latent heat, equal to the different between
the higher internal energy of liquid water and the lower internal energy of solid ice,
at the temperature ``Tᵤ``.

"""
MeltingConstrainedFluxBalance() = MeltingConstrainedFluxBalance(NonlinearSurfaceTemperatureSolver())

#####
##### Flux imbalance and temperature
#####

@inline function top_flux_imbalance(i, j, grid, top_heat_bc, top_surface_temperature,
                                    internal_fluxes, external_fluxes, clock, model_fields)

    # Schematic of the Stefan condition at an upper top
    #
    #  air       ↑   Qx ≡ external_fluxes. Example: Qx = σ T⁴
    #          |⎴⎴⎴|
    # ----------------------- ↕ hₛ(t) → ℒ ∂t hₛ = δQ | T ≤ Tₘ
    #          |⎵⎵⎵|
    #  ice       ↑   Qi ≡ internal_fluxes. Example Qi = - k ∂z T
    #

    Qi = getflux(internal_fluxes, i, j, grid, top_surface_temperature, clock, model_fields)
    Qx = getflux(external_fluxes, i, j, grid, top_surface_temperature, clock, model_fields)

    # The imbalance is defined such that
    # negative imbalance => heat accumulates (out > in) ⟹  growth
    return Qx - Qi
end

@inline top_surface_temperature(i, j, grid, ::PrescribedTemperature, Tu, args...) = Tu

using RootSolvers: SecantMethod, find_zero, CompactSolution

@inline function top_surface_temperature(i, j, grid, top_heat_bc,
                                         current_top_surface_temperature,
                                         internal_fluxes, external_fluxes,
                                         clock, model_fields)

    Tu = get_tracer(i, j, 1, grid, current_top_surface_temperature)
    T₁ = Tu + 1
    T₂ = Tu - 0
    FT = eltype(grid)
    method = SecantMethod{FT}(T₁, T₂)
    solution_type = CompactSolution()

    flux_balance(T) = getflux(external_fluxes, i, j, grid, T, clock, model_fields) -
                      getflux(internal_fluxes, i, j, grid, T, clock, model_fields)

    solution = find_zero(flux_balance, method, solution_type)

    return solution.root
end
