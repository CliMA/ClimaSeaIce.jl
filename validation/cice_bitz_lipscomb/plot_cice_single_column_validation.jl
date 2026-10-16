using CairoMakie

const RESULTS_DIR = joinpath(@__DIR__, "results")
const SUMMARY_CSV = joinpath(RESULTS_DIR, "cice_case_comparison_summary.csv")
const OUTPUT = joinpath(RESULTS_DIR, "cice_single_column_validation.png")

function read_simple_csv(path)
    lines = readlines(path)
    headers = split(lines[1], ','; keepempty = true)
    rows = Vector{Dict{String, String}}()

    for line in lines[2:end]
        isempty(strip(line)) && continue
        values = split(line, ','; keepempty = true)
        row = Dict{String, String}()
        for (header, value) in zip(headers, values)
            row[header] = value
        end
        push!(rows, row)
    end

    return rows
end

function number(row, column)
    value = get(row, column, "")
    return isempty(value) ? NaN : parse(Float64, value)
end

case_order(case_id) =
    case_id == "cold_conductive_relaxation" ? 1 :
    case_id == "surface_warming" ? 2 :
    case_id == "surface_ablation" ? 3 :
    case_id == "basal_growth" ? 4 : 5

case_short_name(case_id) =
    case_id == "cold_conductive_relaxation" ? "cold" :
    case_id == "surface_warming" ? "warming" :
    case_id == "surface_ablation" ? "ablation" :
    case_id == "basal_growth" ? "growth" : case_id

conduct_order(conduct) = conduct == "MU71" ? 1 : 2

rows = sort(read_simple_csv(SUMMARY_CSV);
            by = row -> (case_order(row["case_id"]), conduct_order(row["conduct"])))

labels = [case_short_name(row["case_id"]) * "\n" * row["conduct"] for row in rows]
x = 1:length(rows)

strict_prognostic_acceptance =
    [number(row, "icepack_matrix_prognostic_thickness_forced_replay_acceptance_ratio") for row in rows]

strict_temperature_acceptance =
    [number(row, "icepack_temperature_matrix_cice_thickness_forced_temperature_acceptance_ratio") for row in rows]

strict_enthalpy_percent =
    100 .* [number(row, "icepack_temperature_matrix_cice_thickness_forced_max_relative_cice_enthalpy_error") for row in rows]

strict_column_energy_percent =
    100 .* [number(row, "icepack_temperature_matrix_cice_thickness_forced_max_relative_column_energy_error") for row in rows]

general_temperature_acceptance =
    [number(row, "cice_thickness_forced_temperature_acceptance_ratio") for row in rows]

moving_metric_temperature_acceptance =
    [number(row, "moving_metric_cice_thickness_forced_temperature_acceptance_ratio") for row in rows]

cice_final_hi = [number(row, "final_hi_m") for row in rows]
source_final_hi = [number(row, "icepack_matrix_prognostic_thickness_forced_replay_final_hi_m") for row in rows]

cice_top_melt = [number(row, "integrated_meltt_m") for row in rows]
source_top_melt = [number(row, "icepack_matrix_prognostic_thickness_forced_replay_top_melt_m") for row in rows]

cice_basal_growth = [number(row, "integrated_congel_m") for row in rows]
source_basal_growth = [number(row, "icepack_matrix_prognostic_thickness_forced_replay_basal_growth_m") for row in rows]

colors = Makie.wong_colors()
source_color = colors[1]
temperature_color = colors[2]
general_color = colors[3]
moving_color = colors[4]
enthalpy_color = colors[5]
column_energy_color = colors[6]
final_hi_color = colors[1]
top_melt_color = colors[2]
basal_growth_color = colors[3]

figure = Figure(size = (1500, 980),
                fontsize = 18,
                backgroundcolor = :white)

Label(figure[1, 1:2],
      "ClimaSeaIce single-column BL99 validation against CICE/Icepack",
      fontsize = 26,
      font = :bold,
      tellwidth = false)

Label(figure[2, 1:2],
      "All eight no-snow histories pass the strict source-level BL99 gate. The general and moving-metric replays are shown as diagnostics.",
      fontsize = 16,
      color = (:black, 0.68),
      tellwidth = false)

axis_kwargs = (xticks = (x, labels),
               xticklabelrotation = pi / 6,
               xgridvisible = false)

ax1 = Axis(figure[3, 1];
           title = "Strict pass metrics, acceptance ratio <= 1",
           ylabel = "acceptance ratio",
           axis_kwargs...)

barplot!(ax1, x .- 0.18, strict_prognostic_acceptance;
         width = 0.34, color = source_color, label = "prognostic thickness replay")
barplot!(ax1, x .+ 0.18, strict_temperature_acceptance;
         width = 0.34, color = temperature_color, label = "CICE-thickness temperature")
hlines!(ax1, [1.0]; color = :black, linestyle = :dash, linewidth = 2, label = "1% gate")
ylims!(ax1, 0, 1.15)
axislegend(ax1; position = :lb, framevisible = false)

ax2 = Axis(figure[3, 2];
           title = "Strict source-level energy errors",
           ylabel = "relative error (%)",
           yscale = log10,
           axis_kwargs...)

barplot!(ax2, x .- 0.18, strict_enthalpy_percent;
         width = 0.34, color = enthalpy_color, label = "CICE negative enthalpy")
barplot!(ax2, x .+ 0.18, strict_column_energy_percent;
         width = 0.34, color = column_energy_color, label = "column energy")
hlines!(ax2, [1.0]; color = :black, linestyle = :dash, linewidth = 2, label = "1% gate")
ylims!(ax2, 1e-5, 1.5)
axislegend(ax2; position = :lt, framevisible = false)

ax3 = Axis(figure[4, 1];
           title = "Thickness, melt, and growth components",
           xlabel = "CICE value (m)",
           ylabel = "ClimaSeaIce source-level replay (m)")

scatter!(ax3, cice_final_hi, source_final_hi;
         color = final_hi_color, marker = :circle, markersize = 14, label = "final hi")
scatter!(ax3, cice_top_melt, source_top_melt;
         color = top_melt_color, marker = :rect, markersize = 14, label = "top melt")
scatter!(ax3, cice_basal_growth, source_basal_growth;
         color = basal_growth_color, marker = :utriangle, markersize = 14, label = "basal growth")
lines!(ax3, [0, 1.05], [0, 1.05]; color = :black, linestyle = :dash, linewidth = 2, label = "1:1")
xlims!(ax3, -0.02, 1.08)
ylims!(ax3, -0.02, 1.08)
axislegend(ax3; position = :lt, framevisible = false)

ax4 = Axis(figure[4, 2];
           title = "Diagnostic temperature replays",
           ylabel = "acceptance ratio",
           yscale = log10,
           axis_kwargs...)

barplot!(ax4, x .- 0.24, strict_temperature_acceptance;
         width = 0.22, color = temperature_color, label = "source-level strict")
barplot!(ax4, x, general_temperature_acceptance;
         width = 0.22, color = general_color, label = "general split solve")
barplot!(ax4, x .+ 0.24, moving_metric_temperature_acceptance;
         width = 0.22, color = moving_color, label = "moving metric")
hlines!(ax4, [1.0]; color = :black, linestyle = :dash, linewidth = 2, label = "1% gate")
ylims!(ax4, 1e-4, 1e1)
axislegend(ax4; position = :lt, framevisible = false)

for ax in (ax1, ax2, ax4)
    ax.xticklabelalign = (:right, :center)
end

save(OUTPUT, figure; px_per_unit = 2)
println("Wrote ", OUTPUT)
