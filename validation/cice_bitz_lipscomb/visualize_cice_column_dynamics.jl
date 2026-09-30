const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DOCS_ENV = joinpath(PROJECT_ROOT, "docs")
DOCS_ENV in LOAD_PATH || push!(LOAD_PATH, DOCS_ENV)

using CairoMakie

const VIS_ROOT = @__DIR__
const RESULTS_DIR = joinpath(VIS_ROOT, "results")
const VIS_PREFIX = get(ENV, "CICE_VIS_RESULTS_PREFIX", "basal_growth_MU71")
const VIS_CASE_ID = get(ENV, "CICE_VIS_CASE_ID", "basal_growth")
const VIS_CONDUCT = get(ENV, "CICE_VIS_CONDUCT", "MU71")

haskey(ENV, "CICE_CASE_ID") || (ENV["CICE_CASE_ID"] = VIS_CASE_ID)
haskey(ENV, "CICE_CONDUCT") || (ENV["CICE_CONDUCT"] = VIS_CONDUCT)

include(joinpath(VIS_ROOT, "analyze_cice_comparison.jl"))

const CICE_SUMMARY_PATH = joinpath(RESULTS_DIR, "$(VIS_PREFIX)_cice_smoke_summary.csv")
const FIGURE_PATH = joinpath(RESULTS_DIR, "$(VIS_PREFIX)_column_dynamics.png")
const MOVIE_PATH = joinpath(RESULTS_DIR, "$(VIS_PREFIX)_column_dynamics.mp4")

function read_simple_csv(path)
    lines = readlines(path)
    header = split(first(lines), ','; keepempty = true)
    rows = Vector{Dict{String, String}}()

    for line in lines[2:end]
        isempty(strip(line)) && continue
        values = split(line, ','; keepempty = true)
        push!(rows, Dict(header .=> values))
    end

    return rows
end

parse_float(value) = isempty(value) ? NaN : parse(Float64, value)
parse_int(value) = isempty(value) ? 0 : parse(Int, value)

function parse_profile(row, prefix)
    values = Float64[]
    for k in 1:64
        key = "$(prefix)_$k"
        haskey(row, key) || break
        value = parse_float(row[key])
        isfinite(value) && push!(values, value)
    end
    return values
end

function read_cice_summary_records(path)
    rows = read_simple_csv(path)
    records = CICEHistoryRecord[]

    for row in rows
        push!(records, CICEHistoryRecord(
            source_file = row["source_file"],
            frequency = row["frequency"],
            time_days = parse_float(row["time_days"]),
            active_j = parse_int(row["active_j"]),
            active_i = parse_int(row["active_i"]),
            hi = parse_float(row["hi_m"]),
            hs = parse_float(row["hs_m"]),
            aice = parse_float(row["aice"]),
            Tsfc = parse_float(row["Tsfc_C"]),
            congel = parse_float(row["congel_cm_per_day"]),
            meltt = parse_float(row["meltt_cm_per_day"]),
            meltb = parse_float(row["meltb_cm_per_day"]),
            fsurf_ai = parse_float(row["fsurf_ai_W_m2"]),
            fcondtop_ai = parse_float(row["fcondtop_ai_W_m2"]),
            fbot = parse_float(row["fbot_W_m2"]),
            siflcondbot = parse_float(row["siflcondbot_W_m2"]),
            fhocn = parse_float(row["fhocn_W_m2"]),
            fhocn_ai = parse_float(row["fhocn_ai_W_m2"]),
            fswabs = parse_float(row["fswabs_W_m2"]),
            fswthru = parse_float(row["fswthru_W_m2"]),
            fsens = parse_float(row["fsens_W_m2"]),
            flat = parse_float(row["flat_W_m2"]),
            flwdn = parse_float(row["flwdn_W_m2"]),
            flwup = parse_float(row["flwup_W_m2"]),
            Tair = parse_float(row["Tair_C"]),
            sst = parse_float(row["sst_C"]),
            frzmlt = parse_float(row["frzmlt_W_m2"]),
            Tice = parse_profile(row, "Tice_C"),
            Sice = parse_profile(row, "Sice_ppt"),
            Tsnow = parse_profile(row, "Tsnow_C"),
        ))
    end

    sort!(records; by = record -> record.time_days)
    return records
end

function replay_climaseaice_dynamics(records)
    hourly_records = filter(record -> record.frequency == "hourly", records)
    initial_record = first(records)
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    energy_transport = ConductiveTemperatureTransport(conductivity = cice_conductivity(Float64))
    fixed_thickness = VIS_CASE_ID == "cold_conductive_relaxation"

    current_record = initial_record
    current_temperature = copy(initial_record.Tice)
    current_thickness = initial_record.hi

    times = Float64[initial_record.time_days]
    temperatures = Vector{Vector{Float64}}([copy(initial_record.Tice)])
    thicknesses = Float64[initial_record.hi]

    for target_record in hourly_records
        step = prognostic_thickness_forced_step(current_record,
                                                current_temperature,
                                                current_thickness,
                                                target_record,
                                                relation,
                                                energy_transport;
                                                fixed_thickness)
        current_record = target_record
        current_temperature = copy(step.temperature)
        current_thickness = step.thickness

        push!(times, target_record.time_days)
        push!(temperatures, current_temperature)
        push!(thicknesses, current_thickness)
    end

    return times, temperatures, thicknesses
end

function cice_dynamics(records)
    selected = vcat(first(records), filter(record -> record.frequency == "hourly", records))
    times = [record.time_days for record in selected]
    temperatures = [copy(record.Tice) for record in selected]
    thicknesses = [record.hi for record in selected]
    return times, temperatures, thicknesses
end

function temperature_raster(temperatures, thicknesses, depths)
    image = fill(NaN, length(depths), length(temperatures))

    for n in eachindex(temperatures)
        profile = temperatures[n]
        hi = thicknesses[n]
        layers = length(profile)
        Δz = hi / layers

        for j in eachindex(depths)
            depth = depths[j]
            if depth <= hi + 1e-12
                k = min(layers, floor(Int, depth / Δz) + 1)
                image[j, n] = profile[k]
            end
        end
    end

    return image
end

function step_profile(profile, hi)
    layers = length(profile)
    Δz = hi / layers
    temperatures = Float64[]
    depths = Float64[]

    for k in 1:layers
        top = (k - 1) * Δz
        bottom = k * Δz
        append!(temperatures, (profile[k], profile[k]))
        append!(depths, (top, bottom))
    end

    return temperatures, depths
end

function nearest_index(values, target)
    _, idx = findmin(abs.(values .- target))
    return idx
end

function title_case_name(case_id)
    words = split(replace(case_id, "_" => " "))
    return join(uppercasefirst.(words), ' ')
end

records = read_cice_summary_records(CICE_SUMMARY_PATH)
cice_times, cice_temperatures, cice_thicknesses = cice_dynamics(records)
clima_times, clima_temperatures, clima_thicknesses = replay_climaseaice_dynamics(records)

length(cice_times) == length(clima_times) || error("CICE and ClimaSeaIce trajectories have mismatched lengths.")
all(t -> isapprox(cice_times[t], clima_times[t]; atol = 1e-12, rtol = 0), eachindex(cice_times)) ||
    error("CICE and ClimaSeaIce time grids differ.")

times = cice_times
max_depth = 1.05 * max(maximum(cice_thicknesses), maximum(clima_thicknesses))
depths = collect(range(0, max_depth; length = 240))

cice_image = temperature_raster(cice_temperatures, cice_thicknesses, depths)
clima_image = temperature_raster(clima_temperatures, clima_thicknesses, depths)
ΔT_image = clima_image .- cice_image

finite_temperature_values = filter(isfinite, vcat(vec(cice_image), vec(clima_image)))
temperature_limits = extrema(finite_temperature_values)
finite_delta_values = filter(isfinite, vec(ΔT_image))
delta_limit = max(abs(minimum(finite_delta_values)), abs(maximum(finite_delta_values)))

snapshot_indices = unique([
    1,
    nearest_index(times, 0.5 * last(times)),
    length(times),
])
snapshot_labels = ["day $(round(times[idx]; digits = 1))" for idx in snapshot_indices]
snapshot_colors = [Makie.wong_colors()[1], Makie.wong_colors()[2], Makie.wong_colors()[3]]

figure = Figure(size = (1700, 1150), fontsize = 18, backgroundcolor = :white)

Label(figure[1, 1:2],
      "$(title_case_name(VIS_CASE_ID)) / $(VIS_CONDUCT): CICE vs ClimaSeaIce column dynamics",
      fontsize = 28,
      font = :bold,
      tellwidth = false)

Label(figure[2, 1:2],
      "Hourly layer temperatures are shown against depth below the surface; thickness evolution is plotted separately.",
      fontsize = 16,
      color = (:black, 0.68),
      tellwidth = false)

ax_cice = Axis(figure[3, 1],
               title = "CICE temperature evolution",
               xlabel = "time (days)",
               ylabel = "depth below surface (m)",
               yreversed = true)
hm_cice = heatmap!(ax_cice, times, depths, cice_image;
                   colormap = :thermal,
                   colorrange = temperature_limits,
                   nan_color = :white)
Colorbar(figure[3, 2], hm_cice, label = "temperature (C)")

ax_clima = Axis(figure[4, 1],
                title = "ClimaSeaIce prognostic replay",
                xlabel = "time (days)",
                ylabel = "depth below surface (m)",
                yreversed = true)
heatmap!(ax_clima, times, depths, clima_image;
         colormap = :thermal,
         colorrange = temperature_limits,
         nan_color = :white)

ax_delta = Axis(figure[5, 1],
                title = "ClimaSeaIce minus CICE",
                xlabel = "time (days)",
                ylabel = "depth below surface (m)",
                yreversed = true)
hm_delta = heatmap!(ax_delta, times, depths, ΔT_image;
                    colormap = :balance,
                    colorrange = (-delta_limit, delta_limit),
                    nan_color = :white)
Colorbar(figure[5, 2], hm_delta, label = "delta T (C)")

ax_thickness = Axis(figure[3, 3],
                    title = "Ice thickness evolution",
                    xlabel = "time (days)",
                    ylabel = "thickness (m)")
lines!(ax_thickness, times, cice_thicknesses; linewidth = 4, color = :black, label = "CICE")
lines!(ax_thickness, times, clima_thicknesses; linewidth = 3, color = Makie.wong_colors()[4], label = "ClimaSeaIce")
axislegend(ax_thickness; position = :lt, framevisible = false)

ax_profiles = Axis(figure[4:5, 3],
                   title = "Temperature snapshots",
                   xlabel = "temperature (C)",
                   ylabel = "depth below surface (m)",
                   yreversed = true)

for (snapshot_color, idx, label) in zip(snapshot_colors, snapshot_indices, snapshot_labels)
    cice_profile_temperature, cice_profile_depth = step_profile(cice_temperatures[idx], cice_thicknesses[idx])
    clima_profile_temperature, clima_profile_depth = step_profile(clima_temperatures[idx], clima_thicknesses[idx])
    lines!(ax_profiles, cice_profile_temperature, cice_profile_depth;
           color = snapshot_color, linewidth = 4, label = "CICE $label")
    lines!(ax_profiles, clima_profile_temperature, clima_profile_depth;
           color = snapshot_color, linewidth = 3, linestyle = :dash, label = "ClimaSeaIce $label")
end

axislegend(ax_profiles; position = :rb, framevisible = false)

save(FIGURE_PATH, figure; px_per_unit = 2)

movie = Figure(size = (1700, 1000), fontsize = 18, backgroundcolor = :white)

movie_title = Label(movie[1, 1:2], "", fontsize = 26, font = :bold, tellwidth = false)

movie_ax_cice = Axis(movie[2, 1],
                     title = "CICE",
                     xlabel = "time (days)",
                     ylabel = "depth below surface (m)",
                     yreversed = true)
heatmap!(movie_ax_cice, times, depths, cice_image;
         colormap = :thermal,
         colorrange = temperature_limits,
         nan_color = :white)

movie_ax_clima = Axis(movie[2, 2],
                      title = "ClimaSeaIce",
                      xlabel = "time (days)",
                      ylabel = "depth below surface (m)",
                      yreversed = true)
heatmap!(movie_ax_clima, times, depths, clima_image;
         colormap = :thermal,
         colorrange = temperature_limits,
         nan_color = :white)

movie_ax_thickness = Axis(movie[3, 1],
                          title = "Thickness",
                          xlabel = "time (days)",
                          ylabel = "thickness (m)")
lines!(movie_ax_thickness, times, cice_thicknesses; linewidth = 4, color = :black, label = "CICE")
lines!(movie_ax_thickness, times, clima_thicknesses; linewidth = 3, color = Makie.wong_colors()[4], label = "ClimaSeaIce")
axislegend(movie_ax_thickness; position = :lt, framevisible = false)

movie_ax_profile = Axis(movie[3, 2],
                        title = "Current temperature profile",
                        xlabel = "temperature (C)",
                        ylabel = "depth below surface (m)",
                        yreversed = true)

current_time_line_cice = Observable(Point2f[(times[1], 0), (times[1], max_depth)])
current_time_line_clima = Observable(Point2f[(times[1], 0), (times[1], max_depth)])
lines!(movie_ax_cice, current_time_line_cice; color = :white, linewidth = 3)
lines!(movie_ax_clima, current_time_line_clima; color = :white, linewidth = 3)

current_cice_point = Observable(Point2f[(times[1], cice_thicknesses[1])])
current_clima_point = Observable(Point2f[(times[1], clima_thicknesses[1])])
scatter!(movie_ax_thickness, current_cice_point; color = :black, markersize = 18)
scatter!(movie_ax_thickness, current_clima_point; color = Makie.wong_colors()[4], markersize = 16)

cice_profile_x = Observable(Float64[])
cice_profile_y = Observable(Float64[])
clima_profile_x = Observable(Float64[])
clima_profile_y = Observable(Float64[])
lines!(movie_ax_profile, cice_profile_x, cice_profile_y; color = :black, linewidth = 4, label = "CICE")
lines!(movie_ax_profile, clima_profile_x, clima_profile_y; color = Makie.wong_colors()[4], linewidth = 3, linestyle = :dash, label = "ClimaSeaIce")
axislegend(movie_ax_profile; position = :rb, framevisible = false)

all_temperatures = vcat(reduce(vcat, cice_temperatures), reduce(vcat, clima_temperatures))
x_limits = extrema(all_temperatures)
x_padding = 0.5
xlims!(movie_ax_profile, x_limits[1] - x_padding, x_limits[2] + x_padding)
ylims!(movie_ax_profile, 0, max_depth)
ylims!(movie_ax_cice, 0, max_depth)
ylims!(movie_ax_clima, 0, max_depth)

record(movie, MOVIE_PATH, eachindex(times); framerate = 12) do n
    movie_title.text = "$(title_case_name(VIS_CASE_ID)) / $(VIS_CONDUCT)  day $(round(times[n]; digits = 2))"
    current_time_line_cice[] = Point2f[(times[n], 0), (times[n], max_depth)]
    current_time_line_clima[] = Point2f[(times[n], 0), (times[n], max_depth)]
    current_cice_point[] = Point2f[(times[n], cice_thicknesses[n])]
    current_clima_point[] = Point2f[(times[n], clima_thicknesses[n])]
    cice_x, cice_y = step_profile(cice_temperatures[n], cice_thicknesses[n])
    clima_x, clima_y = step_profile(clima_temperatures[n], clima_thicknesses[n])
    cice_profile_x[] = cice_x
    cice_profile_y[] = cice_y
    clima_profile_x[] = clima_x
    clima_profile_y[] = clima_y
end

println("Wrote ", FIGURE_PATH)
println("Wrote ", MOVIE_PATH)
