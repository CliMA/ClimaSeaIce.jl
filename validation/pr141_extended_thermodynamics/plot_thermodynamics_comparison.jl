#!/usr/bin/env julia

using CairoMakie
using Printf
using TOML

const VALIDATION_DIR = @__DIR__
const DEFAULT_RESULTS_DIR = joinpath(VALIDATION_DIR, "results")
const RESULTS_DIR = get(ENV, "PR141_THERMO_RESULTS_DIR", DEFAULT_RESULTS_DIR)
const SPECS_PATH = joinpath(VALIDATION_DIR, "thermodynamics_case_specs.toml")

function read_csv_rows(path)
    lines = readlines(path)
    header = Symbol.(split(first(lines), ','))
    rows = NamedTuple[]

    for line in lines[2:end]
        isempty(line) && continue
        values = split(line, ',')
        parsed = map(values) do value
            if value == "pass" || value == "fail" || isempty(value)
                return value
            elseif occursin(r"^-?[0-9.]+(?:[eE][-+]?[0-9]+)?$", value)
                return parse(Float64, value)
            else
                return value
            end
        end
        push!(rows, NamedTuple{Tuple(header)}(Tuple(parsed)))
    end

    return rows
end

function spring_spec(specs)
    return only(filter(case -> case["id"] == "melting_in_spring", specs["case"]))
end

function seasonal_spec(specs)
    return only(filter(case -> case["id"] == "arctic_basin_seasonal_cycle", specs["case"]))
end

function plot_melting_in_spring(specs)
    case_spec = spring_spec(specs)
    timeseries_path = joinpath(RESULTS_DIR, case_spec["timeseries_csv"])
    rows = read_csv_rows(timeseries_path)

    forcing_values = sort(unique(row.forcing_wm2 for row in rows))
    colors = Makie.wong_colors()[1:length(forcing_values)]
    labels = [@sprintf("%g", forcing) for forcing in forcing_values]

    fig = Figure(size = (1200, 1000))
    ax_T = Axis(fig[1, 1], ylabel = "Surface temperature (C)",
                title = "Surface temperature: bare (solid) vs snow-covered (dashed)")
    ax_a = Axis(fig[2, 1], ylabel = "Ice concentration (-)",
                title = "Ice concentration: bare (solid) vs snow-covered (dashed)")
    ax_h = Axis(fig[3, 1], ylabel = "Ice thickness (m)",
                title = "Ice thickness: bare (solid) vs snow-covered (dashed)")
    ax_hs = Axis(fig[4, 1], xlabel = "Time (days)", ylabel = "Snow thickness (m)",
                 title = "Snow thickness evolution")

    for (index, forcing) in enumerate(forcing_values)
        bare = sort(filter(row -> row.kind == "bare" && row.forcing_wm2 == forcing, rows); by = row -> row.time_days)
        snow = sort(filter(row -> row.kind == "snow" && row.forcing_wm2 == forcing, rows); by = row -> row.time_days)

        lines!(ax_T, [row.time_days for row in bare], [row.tsfc_c for row in bare], color = colors[index], label = labels[index] * " W/m^2")
        lines!(ax_T, [row.time_days for row in snow], [row.tsfc_c for row in snow], color = colors[index], linestyle = :dash)
        lines!(ax_a, [row.time_days for row in bare], [row.a for row in bare], color = colors[index])
        lines!(ax_a, [row.time_days for row in snow], [row.a for row in snow], color = colors[index], linestyle = :dash)
        lines!(ax_h, [row.time_days for row in bare], [row.h_m for row in bare], color = colors[index])
        lines!(ax_h, [row.time_days for row in snow], [row.h_m for row in snow], color = colors[index], linestyle = :dash)
        lines!(ax_hs, [row.time_days for row in snow], [row.hs_m for row in snow], color = colors[index])
    end

    axislegend(ax_T, position = :rt)
    figure_path = joinpath(RESULTS_DIR, case_spec["figure_png"])
    save(figure_path, fig)
    println("Wrote " * figure_path)
end

function plot_arctic_basin_seasonal_cycle(specs)
    case_spec = seasonal_spec(specs)
    timeseries_path = joinpath(RESULTS_DIR, case_spec["timeseries_csv"])
    rows = sort(read_csv_rows(timeseries_path); by = row -> row.time_days)

    fig = Figure(size = (1000, 1200))
    axT = Axis(fig[1, 1], xlabel = "Time (days)", ylabel = "Top temperature (C)")
    axh = Axis(fig[2, 1], xlabel = "Time (days)", ylabel = "Ice thickness (m)")
    axa = Axis(fig[3, 1], xlabel = "Time (days)", ylabel = "Ice concentration (-)")
    axQ = Axis(fig[4, 1], xlabel = "Time (days)", ylabel = "Top heat flux (W m^-2)")

    lines!(axT, [row.time_days for row in rows], [row.tsfc_c for row in rows])
    lines!(axh, [row.time_days for row in rows], [row.h_m for row in rows])
    lines!(axa, [row.time_days for row in rows], [row.a for row in rows])
    lines!(axQ, [row.time_days for row in rows], [row.q_top_w_m2 for row in rows])

    figure_path = joinpath(RESULTS_DIR, case_spec["figure_png"])
    save(figure_path, fig)
    println("Wrote " * figure_path)
end

function main()
    specs = TOML.parsefile(SPECS_PATH)
    println("PR141 extended thermodynamics plotting")
    println("=====================================")
    println("Spec file: " * SPECS_PATH)
    println("Results dir: " * RESULTS_DIR)
    plot_melting_in_spring(specs)
    plot_arctic_basin_seasonal_cycle(specs)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
