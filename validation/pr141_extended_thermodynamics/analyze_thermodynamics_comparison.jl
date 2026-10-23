#!/usr/bin/env julia

# This validation aims to show that the original example cases `melting_in_spring`
# and `arctic_basin_seasonal_cycle.jl` work well with modest parameter tuning,
# including relative to the existing CICE/BL99 thermodynamic behavior.
# It also writes summary tables and plot-ready timeseries for the standard
# plotting step in this validation workflow.

using Printf
using Statistics
using TOML

using Oceananigans
using Oceananigans.Units
using Oceananigans.Units: Time
using ClimaSeaIce
using ClimaSeaIce.SeaIceThermodynamics.HeatBoundaryConditions: FluxFunction, RadiativeEmission, getflux
const VALIDATION_DIR = @__DIR__
const DEFAULT_RESULTS_DIR = joinpath(VALIDATION_DIR, "results")
const RESULTS_DIR = get(ENV, "PR141_THERMO_RESULTS_DIR", DEFAULT_RESULTS_DIR)
const SPECS_PATH = joinpath(VALIDATION_DIR, "thermodynamics_case_specs.toml")
const DASHBOARD_PATH = joinpath(RESULTS_DIR, "pre_cice_dashboard.csv")
const CASE_COMPARISON_PATH = joinpath(RESULTS_DIR, "pr141_extended_thermodynamics_case_comparison.csv")

mkpath(RESULTS_DIR)

format_bool(x::Bool) = x ? "pass" : "fail"
format_bool(x) = string(x)
format_float(x::AbstractFloat) = @sprintf("%.12g", x)
format_float(x) = string(x)

function csv_field(value)
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    end

    return text
end

function write_csv(path, header, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join(csv_field.(row), ','))
        end
    end
end

function write_named_rows_csv(path, rows)
    field_names = collect(keys(first(rows)))
    table = [[row[name] isa Bool ? format_bool(row[name]) :
              row[name] isa AbstractFloat ? format_float(row[name]) :
              row[name]
              for name in field_names]
             for row in rows]
    write_csv(path, field_names, table)
end

function write_timeseries_csv(path, rows)
    write_named_rows_csv(path, rows)
end

function write_dashboard(path, rows)
    write_csv(path, ("gate", "status", "evidence"),
              [[row.gate, row.status, row.evidence] for row in rows])
end

function first_meltout_day(times, values; threshold = 1e-6)
    for (time_days, value) in zip(times, values)
        value <= threshold && return time_days
    end

    return last(times)
end

function spring_spec(specs)
    return only(filter(case -> case["id"] == "melting_in_spring", specs["case"]))
end

function seasonal_spec(specs)
    return only(filter(case -> case["id"] == "arctic_basin_seasonal_cycle", specs["case"]))
end

function run_melting_in_spring(specs)
    case_spec = spring_spec(specs)
    forcing_values = Float64.(case_spec["forcing_wm2"])
    initial_snow_thickness = Float64(case_spec["initial_snow_thickness_m"])

    grid = RectilinearGrid(size = length(forcing_values), x = (0, 1), topology = (Periodic, Flat, Flat))

    solar_insolation = reshape(forcing_values, (length(forcing_values), 1, 1))
    outgoing_radiation = RadiativeEmission()

    parameters = (
        transfer_coefficient     = 1e-3,
        atmosphere_density       = 1.225,
        atmosphere_heat_capacity = 1004.0,
        atmosphere_temperature   = -5.0,
        atmosphere_wind_speed    = 5.0,
    )

    @inline function sensible_heat_flux(i, j, grid, Tu, clock, fields, parameters)
        Cs = parameters.transfer_coefficient
        rho_a = parameters.atmosphere_density
        c_a = parameters.atmosphere_heat_capacity
        T_a = parameters.atmosphere_temperature
        u_a = parameters.atmosphere_wind_speed
        area_fraction = fields.ℵ[i, j, 1]

        return Cs * rho_a * c_a * u_a * (Tu - T_a) * area_fraction
    end

    aerodynamic_flux = FluxFunction(sensible_heat_flux; parameters)
    top_heat_flux = (outgoing_radiation, solar_insolation, aerodynamic_flux)

    bare_ice_model = SeaIceModel(grid;
                                 ice_consolidation_thickness = 0.05,
                                 top_heat_flux)
    set!(bare_ice_model, h = 1, ℵ = 1)

    snow_thermodynamics = snow_slab_thermodynamics(grid)
    snowy_ice_model = SeaIceModel(grid;
                                  ice_consolidation_thickness = 0.05,
                                  top_heat_flux,
                                  snow_thermodynamics)
    set!(snowy_ice_model, h = 1, ℵ = 1, hs = initial_snow_thickness)

    delta_t = 10minute
    output_interval = 6hours

    bare_series = NamedTuple[]
    snow_series = NamedTuple[]

    function accumulate_bare(sim)
        T = sim.model.ice_thermodynamics.top_surface_temperature
        h = sim.model.ice_thickness
        a = sim.model.ice_concentration

        push!(bare_series, (
            time_days = time(sim) / day,
            h = [h[i, 1, 1] for i in eachindex(forcing_values)],
            a = [a[i, 1, 1] for i in eachindex(forcing_values)],
            T = [T[i, 1, 1] for i in eachindex(forcing_values)],
        ))
    end

    function accumulate_snow(sim)
        model = sim.model
        T = model.snow_thermodynamics.top_surface_temperature
        h = model.ice_thickness
        a = model.ice_concentration
        hs = model.snow_thickness

        push!(snow_series, (
            time_days = time(sim) / day,
            h = [h[i, 1, 1] for i in eachindex(forcing_values)],
            a = [a[i, 1, 1] for i in eachindex(forcing_values)],
            T = [T[i, 1, 1] for i in eachindex(forcing_values)],
            hs = [hs[i, 1, 1] for i in eachindex(forcing_values)],
        ))
    end

    bare_simulation = Simulation(bare_ice_model, Δt = delta_t, stop_time = 30days)
    bare_simulation.callbacks[:save] = Callback(accumulate_bare, TimeInterval(output_interval))
    accumulate_bare(bare_simulation)
    run!(bare_simulation)

    snow_simulation = Simulation(snowy_ice_model, Δt = delta_t, stop_time = 30days)
    snow_simulation.callbacks[:save] = Callback(accumulate_snow, TimeInterval(output_interval))
    accumulate_snow(snow_simulation)
    run!(snow_simulation)

    final_bare = last(bare_series)
    final_snow = last(snow_series)
    summary_rows = NamedTuple[]
    timeseries_rows = NamedTuple[]

    for entry in bare_series
        for column in eachindex(forcing_values)
            push!(timeseries_rows, (
                kind = "bare",
                forcing_wm2 = forcing_values[column],
                time_days = entry.time_days,
                h_m = entry.h[column],
                a = entry.a[column],
                tsfc_c = entry.T[column],
                hs_m = 0.0,
            ))
        end
    end

    for entry in snow_series
        for column in eachindex(forcing_values)
            push!(timeseries_rows, (
                kind = "snow",
                forcing_wm2 = forcing_values[column],
                time_days = entry.time_days,
                h_m = entry.h[column],
                a = entry.a[column],
                tsfc_c = entry.T[column],
                hs_m = entry.hs[column],
            ))
        end
    end

    for column in eachindex(forcing_values)
        t_bare = [entry.time_days for entry in bare_series]
        t_snow = [entry.time_days for entry in snow_series]
        bare_meltout_day = first_meltout_day(t_bare, [entry.h[column] for entry in bare_series])
        snow_meltout_day = first_meltout_day(t_snow, [entry.h[column] for entry in snow_series])

        push!(summary_rows, (
            forcing_wm2 = forcing_values[column],
            bare_final_h_m = final_bare.h[column],
            snow_final_h_m = final_snow.h[column],
            bare_final_a = final_bare.a[column],
            snow_final_a = final_snow.a[column],
            snow_final_hs_m = final_snow.hs[column],
            bare_min_a = minimum(entry.a[column] for entry in bare_series),
            snow_min_a = minimum(entry.a[column] for entry in snow_series),
            bare_meltout_day = bare_meltout_day,
            snow_meltout_day = snow_meltout_day,
            snow_delays_meltout = snow_meltout_day >= bare_meltout_day,
            snow_partially_melts = final_snow.hs[column] < initial_snow_thickness,
            snow_preserves_more_concentration = final_snow.a[column] >= final_bare.a[column],
        ))
    end

    summary_path = joinpath(RESULTS_DIR, case_spec["summary_csv"])
    write_named_rows_csv(summary_path, summary_rows)
    timeseries_path = joinpath(RESULTS_DIR, case_spec["timeseries_csv"])
    write_timeseries_csv(timeseries_path, timeseries_rows)

    all_delay = all(row.snow_delays_meltout for row in summary_rows)
    all_partial_melt = all(row.snow_partially_melts for row in summary_rows)
    all_concentration = all(row.snow_preserves_more_concentration for row in summary_rows)

    case_row = (
        case_id = case_spec["id"],
        description = case_spec["description"],
        status = all((all_delay, all_partial_melt, all_concentration)) ? "pass" : "fail",
        mean_meltout_delay_days = mean(row.snow_meltout_day - row.bare_meltout_day for row in summary_rows),
        minimum_meltout_delay_days = minimum(row.snow_meltout_day - row.bare_meltout_day for row in summary_rows),
        snow_delays_meltout_status = format_bool(all_delay),
        snow_partially_melts_status = format_bool(all_partial_melt),
        snow_preserves_more_concentration_status = format_bool(all_concentration),
        summary_csv = summary_path,
        timeseries_csv = timeseries_path,
        figure_png = joinpath(RESULTS_DIR, case_spec["figure_png"]),
    )

    return case_row
end

function run_arctic_basin_seasonal_cycle(specs)
    case_spec = seasonal_spec(specs)

    grid = RectilinearGrid(size = (), topology = (Flat, Flat, Flat))

    tabulated_shortwave = -[0, 0, 1.9, 9.9, 17.7, 19.2, 13.6, 9.0, 3.7, 0.4, 0, 0] .* 1e4
    tabulated_longwave  = -[10.4, 10.3, 10.3, 11.6, 15.1, 18.0, 19.1, 18.7, 16.5, 13.9, 11.2, 10.9] .* 1e4
    tabulated_sensible  = -[1.18, 0.76, 0.72, 0.29, -0.45, -0.39, -0.30, -0.40, -0.17, 0.1, 0.56, 0.79] .* 1e4
    tabulated_latent    = -[0, -0.02, -0.03, -0.09, -0.46, -0.70, -0.64, -0.66, -0.39, -0.19, -0.01, -0.01] .* 1e4

    nmonths = 12
    month_days = 30
    year_days = month_days * nmonths
    times_days = collect(15:30:(year_days - 15))
    times = times_days .* day

    emissivity = 1.0
    kcal_to_joules = 4184.0

    tabulated_shortwave .*= kcal_to_joules / (month_days * days)
    tabulated_longwave  .*= kcal_to_joules / (month_days * days) .* emissivity
    tabulated_sensible  .*= kcal_to_joules / (month_days * days)
    tabulated_latent    .*= kcal_to_joules / (month_days * days)

    Rs = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Oceananigans.OutputReaders.Cyclical())
    Rl = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Oceananigans.OutputReaders.Cyclical())
    Qs = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Oceananigans.OutputReaders.Cyclical())
    Ql = FieldTimeSeries{Nothing, Nothing, Nothing}(grid, times; time_indexing = Oceananigans.OutputReaders.Cyclical())

    for (i, _) in enumerate(times)
        set!(Rs[i], tabulated_shortwave[i:i])
        set!(Rl[i], tabulated_longwave[i:i])
        set!(Qs[i], tabulated_sensible[i:i])
        set!(Ql[i], tabulated_latent[i:i])
    end

    @inline function linearly_interpolate_flux(i, j, grid, Ts, clock, model_fields, flux)
        t = Time(clock.time)
        return flux[i, j, 1, t]
    end

    @inline function linearly_interpolate_solar_flux(i, j, grid, Ts, clock, model_fields, flux)
        Q = linearly_interpolate_flux(i, j, grid, Ts, clock, model_fields, flux)
        albedo = ifelse(Ts < -0.1, 0.75, 0.64)
        return Q * (1 - albedo)
    end

    sigma = 5.67e-8 * 1.02

    Q_shortwave = FluxFunction(linearly_interpolate_solar_flux, parameters = Rs)
    Q_longwave  = FluxFunction(linearly_interpolate_flux, parameters = Rl)
    Q_sensible  = FluxFunction(linearly_interpolate_flux, parameters = Qs)
    Q_latent    = FluxFunction(linearly_interpolate_flux, parameters = Ql)
    Q_emission  = RadiativeEmission(emissivity = emissivity, stefan_boltzmann_constant = sigma)
    top_heat_flux = (Q_shortwave, Q_longwave, Q_sensible, Q_latent, Q_emission)

    model = SeaIceModel(grid; top_heat_flux)
    set!(model, h = 0.3, ℵ = 1)

    years_simulated = Int(case_spec["years_simulated"])
    simulation = Simulation(model, Δt = 8hours, stop_time = years_simulated * 360days)

    series = NamedTuple[]
    function accumulate_timeseries(sim)
        T = model.ice_thermodynamics.top_surface_temperature
        h = model.ice_thickness
        a = model.ice_concentration
        Qe = model.external_heat_fluxes.top
        Qe = getflux(Qe, 1, 1, grid, first(T), model.clock, model.ice_thickness)

        push!(series, (
            time_days = time(sim) / day,
            h = first(h),
            T = first(T),
            a = first(a),
            Q = Qe,
        ))
    end

    simulation.callbacks[:save] = Callback(accumulate_timeseries, TimeInterval(1days))
    accumulate_timeseries(simulation)
    run!(simulation)

    t = [entry.time_days for entry in series]

    final_year_start = last(t) - year_days
    previous_year_start = final_year_start - year_days

    final_year = [entry for entry in series if final_year_start <= entry.time_days <= last(t)]
    previous_year = [entry for entry in series if previous_year_start <= entry.time_days < final_year_start]

    sample_count = min(length(final_year), length(previous_year))
    final_h = [final_year[i].h for i in 1:sample_count]
    previous_h = [previous_year[end - sample_count + i].h for i in 1:sample_count]
    final_T = [final_year[i].T for i in 1:sample_count]
    previous_T = [previous_year[end - sample_count + i].T for i in 1:sample_count]

    max_year_to_year_h_diff = maximum(abs.(final_h .- previous_h))
    max_year_to_year_T_diff = maximum(abs.(final_T .- previous_T))
    thickness_range = maximum(final_h) - minimum(final_h)

    summary_row = (
        years_simulated = years_simulated,
        final_year_h_min_m = minimum(final_h),
        final_year_h_max_m = maximum(final_h),
        final_year_a_min = minimum(entry.a for entry in final_year),
        final_year_a_max = maximum(entry.a for entry in final_year),
        final_year_T_min_C = minimum(final_T),
        final_year_T_max_C = maximum(final_T),
        final_year_Q_min_W_m2 = minimum(entry.Q for entry in final_year),
        final_year_Q_max_W_m2 = maximum(entry.Q for entry in final_year),
        max_year_to_year_h_diff_m = max_year_to_year_h_diff,
        max_year_to_year_T_diff_C = max_year_to_year_T_diff,
        seasonal_cycle_present = thickness_range > Float64(case_spec["minimum_seasonal_thickness_range_m"]),
        quasi_equilibrated_thickness = max_year_to_year_h_diff < Float64(case_spec["maximum_year_to_year_thickness_diff_m"]),
        quasi_equilibrated_temperature = max_year_to_year_T_diff < Float64(case_spec["maximum_year_to_year_temperature_diff_C"]),
    )

    summary_path = joinpath(RESULTS_DIR, case_spec["summary_csv"])
    write_named_rows_csv(summary_path, [summary_row])
    timeseries_path = joinpath(RESULTS_DIR, case_spec["timeseries_csv"])
    write_timeseries_csv(timeseries_path, [(
        time_days = entry.time_days,
        h_m = entry.h,
        a = entry.a,
        tsfc_c = entry.T,
        q_top_w_m2 = entry.Q,
    ) for entry in series])

    case_row = (
        case_id = case_spec["id"],
        description = case_spec["description"],
        status = all((summary_row.seasonal_cycle_present,
                      summary_row.quasi_equilibrated_thickness,
                      summary_row.quasi_equilibrated_temperature)) ? "pass" : "fail",
        final_year_thickness_range_m = thickness_range,
        max_year_to_year_h_diff_m = max_year_to_year_h_diff,
        max_year_to_year_T_diff_C = max_year_to_year_T_diff,
        seasonal_cycle_present_status = format_bool(summary_row.seasonal_cycle_present),
        quasi_equilibrated_thickness_status = format_bool(summary_row.quasi_equilibrated_thickness),
        quasi_equilibrated_temperature_status = format_bool(summary_row.quasi_equilibrated_temperature),
        summary_csv = summary_path,
        timeseries_csv = timeseries_path,
        figure_png = joinpath(RESULTS_DIR, case_spec["figure_png"]),
    )

    return case_row
end

function write_case_comparison(path, case_rows)
    header = Symbol[]

    for row in case_rows
        for name in keys(row)
            name in header || push!(header, name)
        end
    end

    table = [[get(row, name, "") for name in header] for row in case_rows]
    write_csv(path, header, table)
end

function print_case_status(case_row)
    println(case_row.case_id * "_status=" * case_row.status)
    println("  summary_csv=" * case_row.summary_csv)
    println("  figure_png=" * case_row.figure_png)
end

function main()
    specs = TOML.parsefile(SPECS_PATH)

    println("PR141 extended thermodynamics validation")
    println("=======================================")
    println("Spec file: " * SPECS_PATH)
    println("Results dir: " * RESULTS_DIR)

    spring_case = run_melting_in_spring(specs)
    seasonal_case = run_arctic_basin_seasonal_cycle(specs)
    case_rows = [spring_case, seasonal_case]

    write_case_comparison(CASE_COMPARISON_PATH, case_rows)
    write_dashboard(DASHBOARD_PATH,
                    [(gate = row.case_id, status = row.status, evidence = row.summary_csv) for row in case_rows])

    println()
    println("Validation gates")
    println("================")
    for row in case_rows
        @printf("%-32s %-8s %s\n", row.case_id, row.status, row.summary_csv)
    end

    println()
    println("Aggregate comparison: " * CASE_COMPARISON_PATH)
    println("Dashboard: " * DASHBOARD_PATH)
    println()
    print_case_status(spring_case)
    print_case_status(seasonal_case)

    return all(row.status == "pass" for row in case_rows) ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
