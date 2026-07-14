using ClimaSeaIce
using Oceananigans
using Test

# The thermodynamic kernels overwrite `model.mass_fluxes` (kg m⁻² s⁻¹) on every call.
# Mass closure identity (thermodynamics only, no advection):
#
#     ∂t(ρᵢ h ℵ + ρₛ hs ℵ) = thermodynamics.ice + thermodynamics.snow + intercepted_snowfall
#
# With the SplitRungeKutta3 stepper every substep re-bases the prognostic state from Ψ⁻ and the final
# substep advances the full Δt, so the fluxes recorded by the last call match the full-step mass change.

column_value(f) = f[1, 1, 1]

function column_masses(model)
    ρi = column_value(model.sea_ice_density)
    h  = column_value(model.ice_thickness)
    ℵ  = column_value(model.ice_concentration)
    Mi = ρi * h * ℵ

    Ms = if isnothing(model.snow_thickness)
        zero(Mi)
    else
        ρs = column_value(model.snow_density)
        hs = column_value(model.snow_thickness)
        ρs * hs * ℵ
    end

    return Mi, Ms
end

function mass_flux_closure(model, Δt)
    Mi⁻, Ms⁻ = column_masses(model)
    time_step!(model, Δt)
    Mi⁺, Ms⁺ = column_masses(model)

    fluxes = model.mass_fluxes
    total = column_value(fluxes.thermodynamics.ice) +
            column_value(fluxes.thermodynamics.snow) +
            column_value(fluxes.intercepted_snowfall)
    expected = ((Mi⁺ + Ms⁺) - (Mi⁻ + Ms⁻)) / Δt

    return total, expected
end

# Absolute tolerance for the closure identity, scaled by the flux magnitude. The 1e-12 floor assumes
# Float64 grids; for lower-precision grids it should scale with `eps(eltype(grid))`.
closure_tolerance(expected) = 1e-12 * max(1, abs(expected))

@testset "Thermodynamic mass fluxes: bare ice [$timestepper]" for timestepper in (:ForwardEuler, :SplitRungeKutta3)
    Δt = 3600

    @testset "Freezing" begin
        grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
        model = SeaIceModel(grid; top_heat_flux = 100, bottom_heat_flux = 10, timestepper)
        set!(model, h=1, ℵ=1)

        total, expected = mass_flux_closure(model, Δt)
        fluxes = model.mass_fluxes

        @test total ≈ expected atol=closure_tolerance(expected)
        @test column_value(fluxes.thermodynamics.ice) > 0
        @test column_value(fluxes.thermodynamics.snow) == 0
        @test column_value(fluxes.intercepted_snowfall) == 0
    end

    @testset "Melting" begin
        grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
        model = SeaIceModel(grid; top_heat_flux = -200, bottom_heat_flux = 10, timestepper)
        set!(model, h=1, ℵ=1)

        total, expected = mass_flux_closure(model, Δt)
        fluxes = model.mass_fluxes

        @test total ≈ expected atol=closure_tolerance(expected)
        @test column_value(fluxes.thermodynamics.ice) < 0
        @test column_value(fluxes.thermodynamics.snow) == 0
        @test column_value(fluxes.intercepted_snowfall) == 0
    end

    @testset "Melt to extinction" begin
        grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
        model = SeaIceModel(grid; top_heat_flux = -1e5, bottom_heat_flux = 10, timestepper)
        set!(model, h=0.2, ℵ=1)

        total, expected = mass_flux_closure(model, Δt)

        @test column_value(model.ice_thickness) == 0
        @test column_value(model.ice_concentration) == 0
        @test total ≈ expected atol=closure_tolerance(expected)
    end

    @testset "Partial concentration freezing" begin
        # Starting below full cover exercises the lateral-growth concentration path (ℵ increases).
        grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
        model = SeaIceModel(grid; top_heat_flux = 300, bottom_heat_flux = 10, timestepper)
        set!(model, h=1, ℵ=0.95)

        total, expected = mass_flux_closure(model, Δt)
        fluxes = model.mass_fluxes

        @test column_value(model.ice_concentration) > 0.95
        @test column_value(fluxes.thermodynamics.ice) > 0
        @test total ≈ expected atol=closure_tolerance(expected)
    end
end

@testset "Thermodynamic mass fluxes: snow + ice [$timestepper]" for timestepper in (:ForwardEuler, :SplitRungeKutta3)
    Δt = 3600
    Ps = 1e-5 # kg m⁻² s⁻¹ snowfall mass flux

    @testset "Freezing" begin
        grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
        model = SeaIceModel(grid;
                            snow_thermodynamics = snow_slab_thermodynamics(grid),
                            snowfall = Ps,
                            top_heat_flux = 100,
                            bottom_heat_flux = 10,
                            timestepper)
        set!(model, h=1, ℵ=1, hs=0.1)

        total, expected = mass_flux_closure(model, Δt)
        fluxes = model.mass_fluxes
        ℵ⁺ = column_value(model.ice_concentration)

        @test total ≈ expected atol=closure_tolerance(expected)
        @test column_value(fluxes.intercepted_snowfall) ≈ Ps * ℵ⁺ atol=1e-12 * Ps
        @test column_value(fluxes.thermodynamics.ice) > 0
    end

    @testset "Melting" begin
        grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
        model = SeaIceModel(grid;
                            snow_thermodynamics = snow_slab_thermodynamics(grid),
                            snowfall = Ps,
                            top_heat_flux = -200,
                            bottom_heat_flux = 10,
                            timestepper)
        set!(model, h=1, ℵ=1, hs=0.1)

        total, expected = mass_flux_closure(model, Δt)
        fluxes = model.mass_fluxes
        ℵ⁺ = column_value(model.ice_concentration)

        @test total ≈ expected atol=closure_tolerance(expected)
        @test column_value(fluxes.intercepted_snowfall) ≈ Ps * ℵ⁺ atol=1e-12 * Ps
        @test column_value(fluxes.thermodynamics.snow) < 0 # top melt removes snow faster than snowfall replenishes it
    end

    @testset "Melt to extinction" begin
        grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
        model = SeaIceModel(grid;
                            snow_thermodynamics = snow_slab_thermodynamics(grid),
                            snowfall = Ps,
                            top_heat_flux = -1e5,
                            bottom_heat_flux = 10,
                            timestepper)
        set!(model, h=0.2, ℵ=1, hs=0.05)

        total, expected = mass_flux_closure(model, Δt)

        @test column_value(model.ice_thickness) == 0
        @test column_value(model.ice_concentration) == 0
        @test column_value(model.snow_thickness) == 0
        @test total ≈ expected atol=closure_tolerance(expected)
    end
end

# The thermodynamic kernels write the diagnostics on every column, including immersed (land) cells;
# `update_state!` must mask them so land reports zero exchange rather than a phantom flux.
@testset "Thermodynamic mass fluxes: immersed masking [$timestepper]" for timestepper in (:ForwardEuler, :SplitRungeKutta3)
    Δt = 3600

    underlying_grid = RectilinearGrid(size=(2, 1, 1), x=(0, 2), y=(0, 1), z=(-1, 0),
                                      topology=(Bounded, Bounded, Bounded))

    # Column i=1 is land (bottom at the surface), column i=2 is ocean.
    bottom(x, y) = x < 1 ? 0.0 : -1.0
    grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom))

    model = SeaIceModel(grid; top_heat_flux = 100, bottom_heat_flux = 10, timestepper)
    set!(model, h=1, ℵ=1)

    time_step!(model, Δt)

    fluxes = model.mass_fluxes
    ice  = interior(fluxes.thermodynamics.ice)
    snow = interior(fluxes.thermodynamics.snow)
    snowfall = interior(fluxes.intercepted_snowfall)

    @test ice[1, 1, 1]  == 0   # land column masked
    @test snow[1, 1, 1] == 0
    @test snowfall[1, 1, 1] == 0
    @test ice[2, 1, 1]  != 0   # ocean column freezes
end
