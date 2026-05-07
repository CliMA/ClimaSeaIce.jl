using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics: latent_heat, bottom_temperature
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

    snow_thermo = snow ? snow_slab_thermodynamics(grid) : nothing
    Ps = precipitation ? 6e-5 : 0

    model = SeaIceModel(grid;
                        ice_consolidation_thickness = 0.05,
                        top_heat_flux,
                        bottom_heat_flux = bot_heat_flux,
                        snow_thermodynamics = snow_thermo,
                        snowfall = Ps)

    if snow
        set!(model, h=1.0, ℵ=1, hs=0.2)
    else
        set!(model, h=1.0, ℵ=1)
    end

    pt = model.phase_transitions
    ℒ  = latent_heat(pt, 0)                    # per-mass at 0 ᵒC
    ρi = @inbounds model.sea_ice_density[1, 1, 1]  # bulk sea-ice density
    ρs = snow ? @inbounds(model.snow_density[1, 1, 1]) : 0.0

    Δt = 600.0
    Nsteps = 200
    max_residual = 0.0

    for n in 1:Nsteps
        h₀  = first(interior(model.ice_thickness))
        ℵ₀  = first(interior(model.ice_concentration))
        hs₀ = snow ? first(interior(model.snow_thickness)) : 0.0
        E₀  = -ℵ₀ * (ρi * ℒ * h₀ + ρs * ℒ * hs₀)

        time_step!(model, Δt)

        h₁  = first(interior(model.ice_thickness))
        ℵ₁  = first(interior(model.ice_concentration))
        hs₁ = snow ? first(interior(model.snow_thickness)) : 0.0
        E₁  = -ℵ₁ * (ρi * ℒ * h₁ + ρs * ℒ * hs₁)

        Qa = top_record[1]
        Ql = bot_record[1]
        Qp = (precipitation && ℵ₁ > 0) ? -ℒ * Ps : 0.0

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

# The tests above run with ℵ = 1 throughout, which hides a per-cell /
# per-ice convention inconsistency in the layered snow-ice kernel (the
# snow-surface energy balance compared per-cell atmospheric flux with the
# per-ice conductive flux, under-estimating snow melt when ℵ < 1). The
# tests below exercise ℵ < 1 to catch regressions on that fix. They use
# a per-cell top heat flux (consistent with what `_layered_thermodynamic_time_step!`
# actually expects after the upstream NumericalEarth assembler fix).

@inline function recording_per_cell_flux(i, j, grid, Tu, clock, fields, p)
    ℵ = fields.ℵ[i, j, 1]
    # Per-ice flux times ℵ → per-cell (mirrors what the coupled assembler writes)
    Q_per_ice = p.coefficient * (Tu - p.temperature)
    Q         = Q_per_ice * ℵ
    p.record[1] = Q
    return Q
end

function partial_cover_energy_conservation_test(; ℵ₀ = 0.5, hs₀ = 0.15, melting = true)
    grid = RectilinearGrid(size=(), topology=(Flat, Flat, Flat))

    Ta = melting ? 5.0 : -15.0
    top_record = [0.0]
    top_params = (coefficient = 1e-3 * 1.225 * 1004 * 5, temperature = Ta, record = top_record)
    top_heat_flux = FluxFunction(recording_per_cell_flux; parameters = top_params)

    Qb = melting ? -20.0 : -5.0
    bot_record = [0.0]
    bot_heat_flux = FluxFunction(bottom_recording_flux; parameters = (flux_value = Qb, record = bot_record))

    snow_thermo = snow_slab_thermodynamics(grid)

    model = SeaIceModel(grid;
                        ice_consolidation_thickness = 0.05,
                        top_heat_flux,
                        bottom_heat_flux = bot_heat_flux,
                        snow_thermodynamics = snow_thermo,
                        snowfall = 0)

    set!(model, h = 1.0, ℵ = ℵ₀, hs = hs₀)

    pt = model.phase_transitions
    ℒ  = latent_heat(pt, 0)
    ρi = @inbounds model.sea_ice_density[1, 1, 1]
    ρs = @inbounds model.snow_density[1, 1, 1]

    Δt = 600.0
    Nsteps = 200
    max_residual = 0.0

    for n in 1:Nsteps
        h₀  = first(interior(model.ice_thickness))
        ℵ₀_ = first(interior(model.ice_concentration))
        hs₀_ = first(interior(model.snow_thickness))
        E₀  = -ℵ₀_ * (ρi * ℒ * h₀ + ρs * ℒ * hs₀_)

        time_step!(model, Δt)

        h₁  = first(interior(model.ice_thickness))
        ℵ₁  = first(interior(model.ice_concentration))
        hs₁ = first(interior(model.snow_thickness))
        E₁  = -ℵ₁ * (ρi * ℒ * h₁ + ρs * ℒ * hs₁)

        Qa = top_record[1]            # per-cell (mirrors coupler output)
        Ql = bot_record[1]            # per-cell bottom

        dE       = E₁ - E₀
        expected = (-Qa + Ql) * Δt
        scale    = max(abs(E₀), abs(E₁), abs(expected), 1.0)
        residual = abs(dE - expected) / scale
        max_residual = max(max_residual, residual)

        h₁ ≤ 0 && ℵ₁ ≤ 0 && break
    end

    return max_residual
end

@testset "Energy conservation with partial ice cover (ℵ < 1)" begin
    # The closed-form solve in `_layered_thermodynamic_time_step!` is
    # self-consistent to floating-point precision, so the per-step
    # energy residual is at round-off. Use a conservative rtol that
    # still catches the O(1-ℵ) per-cell/per-ice bug (~1e-4 before the
    # fix) but allows room for floating-point accumulation over many
    # steps.
    rtol = 1e-13

    @testset "ℵ = 0.5, snow-covered melting" begin
        @test partial_cover_energy_conservation_test(ℵ₀=0.5, hs₀=0.15, melting=true) < rtol
    end

    @testset "ℵ = 0.8, snow-covered melting" begin
        @test partial_cover_energy_conservation_test(ℵ₀=0.8, hs₀=0.15, melting=true) < rtol
    end

    @testset "ℵ = 0.3, snow-covered freezing" begin
        @test partial_cover_energy_conservation_test(ℵ₀=0.3, hs₀=0.05, melting=false) < rtol
    end
end
