using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: latent_heat
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: FluxFunction
using Oceananigans
using Oceananigans.Fields: interior
using Test

@inline function recording_flux(i, j, grid, Tu, clock, fields, p)
    ℵ = fields.ℵ[i, j, 1]
    Q = p.coefficient * (Tu - p.temperature) * ℵ
    p.record[1] = Q
    return Q
end

@inline function bottom_recording_flux(i, j, grid, Tu, clock, fields, p)
    p.record[1] = p.flux_value
    return p.flux_value
end

function energy_conservation_test(; snow=false, precipitation=false, melting=false)
    grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

    Ta = melting ? 5.0 : -15.0
    top_record = [0.0]
    top_params = (coefficient = 1e-3 * 1.225 * 1004 * 5, temperature = Ta, record = top_record)
    top_heat_flux = FluxFunction(recording_flux; parameters = top_params)

    Qb = melting ? -20.0 : -5.0
    bot_record = [0.0]
    bot_heat_flux = FluxFunction(bottom_recording_flux; parameters = (flux_value = Qb, record = bot_record))

    snow_thermo = snow ? SlabSnowThermodynamics(grid) : nothing
    Ps = precipitation ? 6e-5 : 0

    model = SeaIceModel(grid;
                        ice_consolidation_thickness = 0.05,
                        top_heat_flux,
                        bottom_heat_flux = bot_heat_flux,
                        snow_thermodynamics = snow_thermo,
                        snow_precipitation = Ps)

    if snow
        set!(model, h=1.0, ℵ=1, hs=0.2)
    else
        set!(model, h=1.0, ℵ=1)
    end

    PT = model.ice_thermodynamics.phase_transitions
    ℰ  = latent_heat(PT, 0)
    ρs = snow ? snow_thermo.phase_transitions.density : 0.0
    Ls = snow ? snow_thermo.phase_transitions.reference_latent_heat : 0.0

    Δt = 600.0
    Nsteps = 200
    max_residual = 0.0

    for n in 1:Nsteps
        h₀  = first(interior(model.ice_thickness))
        ℵ₀  = first(interior(model.ice_concentration))
        hs₀ = snow ? first(interior(model.snow_thickness)) : 0.0
        E₀  = -ℵ₀ * (ℰ * h₀ + ρs * Ls * hs₀)

        time_step!(model, Δt)

        h₁  = first(interior(model.ice_thickness))
        ℵ₁  = first(interior(model.ice_concentration))
        hs₁ = snow ? first(interior(model.snow_thickness)) : 0.0
        E₁  = -ℵ₁ * (ℰ * h₁ + ρs * Ls * hs₁)

        Qa = top_record[1]
        Ql = bot_record[1]
        Qp = (precipitation && ℵ₁ > 0) ? -Ls * Ps : 0.0

        dE = E₁ - E₀
        expected = (-Qa + Ql + Qp) * Δt
        scale = max(abs(E₀), abs(E₁), abs(expected), 1.0)
        residual = abs(dE - expected) / scale
        max_residual = max(max_residual, residual)

        h₁ ≤ 0 && ℵ₁ ≤ 0 && break
    end

    return max_residual
end

@testset "Energy conservation" begin
    rtol = 1e-15

    @testset "Bare ice, freezing" begin
        @test energy_conservation_test(snow=false, melting=false) < rtol
    end

    @testset "Bare ice, melting" begin
        @test energy_conservation_test(snow=false, melting=true) < rtol
    end

    @testset "Snow-covered, freezing" begin
        @test energy_conservation_test(snow=true, melting=false) < rtol
    end

    @testset "Snow-covered, melting" begin
        @test energy_conservation_test(snow=true, melting=true) < rtol
    end

    @testset "Snow with precipitation, freezing" begin
        @test energy_conservation_test(snow=true, precipitation=true, melting=false) < rtol
    end

    @testset "Snow with precipitation, melting" begin
        @test energy_conservation_test(snow=true, precipitation=true, melting=true) < rtol
    end
end
