using ClimaSeaIce
using Oceananigans
using Test

# The thermodynamic kernels overwrite `model.thermodynamic_mass_fluxes` (kg m⁻² s⁻¹) on every call.
# Mass closure identity (thermodynamics only, no advection):
#
#     ∂t(ρᵢ h ℵ + ρₛ hs ℵ) = ice + snow + intercepted_snowfall
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

    fluxes = model.thermodynamic_mass_fluxes
    total = column_value(fluxes.ice) + column_value(fluxes.snow) + column_value(fluxes.intercepted_snowfall)
    expected = ((Mi⁺ + Ms⁺) - (Mi⁻ + Ms⁻)) / Δt

    return total, expected
end

closure_tolerance(expected) = 1e-12 * max(1, abs(expected))

@testset "Thermodynamic mass fluxes: bare ice [$timestepper]" for timestepper in (:ForwardEuler, :SplitRungeKutta3)
    Δt = 3600

    @testset "Freezing" begin
        grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
        model = SeaIceModel(grid; top_heat_flux = 100, bottom_heat_flux = 10, timestepper)
        set!(model, h=1, ℵ=1)

        total, expected = mass_flux_closure(model, Δt)
        fluxes = model.thermodynamic_mass_fluxes

        @test total ≈ expected atol=closure_tolerance(expected)
        @test column_value(fluxes.ice) > 0
        @test column_value(fluxes.snow) == 0
        @test column_value(fluxes.intercepted_snowfall) == 0
    end

    @testset "Melting" begin
        grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))
        model = SeaIceModel(grid; top_heat_flux = -200, bottom_heat_flux = 10, timestepper)
        set!(model, h=1, ℵ=1)

        total, expected = mass_flux_closure(model, Δt)
        fluxes = model.thermodynamic_mass_fluxes

        @test total ≈ expected atol=closure_tolerance(expected)
        @test column_value(fluxes.ice) < 0
        @test column_value(fluxes.snow) == 0
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
        fluxes = model.thermodynamic_mass_fluxes
        ℵ⁺ = column_value(model.ice_concentration)

        @test total ≈ expected atol=closure_tolerance(expected)
        @test column_value(fluxes.intercepted_snowfall) ≈ Ps * ℵ⁺ atol=1e-12 * Ps
        @test column_value(fluxes.ice) > 0
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
        fluxes = model.thermodynamic_mass_fluxes
        ℵ⁺ = column_value(model.ice_concentration)

        @test total ≈ expected atol=closure_tolerance(expected)
        @test column_value(fluxes.intercepted_snowfall) ≈ Ps * ℵ⁺ atol=1e-12 * Ps
        @test column_value(fluxes.snow) < 0 # top melt removes snow faster than snowfall replenishes it
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
