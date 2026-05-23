#!/usr/bin/env julia

using Printf
using Oceananigans
using Oceananigans.Fields: interior, set!
using Oceananigans.Grids: MutableVerticalDiscretization
using ClimaSeaIce.SeaIceThermodynamics:
    BubblyBrineConductivity,
    column_energy_thickness_remap!,
    ColumnBoundaryConditions,
    ConductiveTemperatureTransport,
    FixedDrainedIceSalinityProfile,
    FixedSalinityBrinePocketEnergyRelation,
    MaykutUntersteinerConductivity,
    PrescribedEnergyFlux,
    PrescribedTemperature,
    column_energy_time_step!,
    compute_column_thermodynamic_diagnostics!,
    compute_column_transport_coefficients!,
    icepack_temperature_matrix_step!,
    internal_energy,
    ice_thermal_conductivity,
    prescribed_salinity_enthalpy_thermodynamics,
    salinity_at_normalized_depth,
    temperature

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS_DIR = joinpath(@__DIR__, "results")
const DASHBOARD = joinpath(RESULTS_DIR, "pre_cice_dashboard.csv")
const RESULT_PREFIX = replace(get(ENV, "CICE_RESULTS_PREFIX", ""), r"[^A-Za-z0-9_.-]" => "_")
const CICE_SUMMARY = joinpath(RESULTS_DIR, RESULT_PREFIX == "" ? "cice_smoke_summary.csv" :
                                                   "$(RESULT_PREFIX)_cice_smoke_summary.csv")
const CLIMASEAICE_REPLAY_SUMMARY = joinpath(RESULTS_DIR, RESULT_PREFIX == "" ? "climaseaice_state_replay_summary.csv" :
                                                                 "$(RESULT_PREFIX)_climaseaice_state_replay_summary.csv")
const FORCED_REPLAY_SUMMARY = joinpath(RESULTS_DIR, RESULT_PREFIX == "" ? "fixed_grid_forced_replay_summary.csv" :
                                                            "$(RESULT_PREFIX)_fixed_grid_forced_replay_summary.csv")
const CASE_COMPARISON_SUMMARY = joinpath(RESULTS_DIR, "cice_case_comparison_summary.csv")
const DEFAULT_CICE_HISTORY_DIR =
    "/private/tmp/cice-validation/runs/bl99_cold_conductive_relaxation_MU71/history"
const CASE_ID = get(ENV, "CICE_CASE_ID", "")
const CONDUCTIVITY_VARIANT = get(ENV, "CICE_CONDUCT", "MU71")
const INTERNAL_ENERGY_RELATIVE_FLOOR =
    parse(Float64, get(ENV, "CICE_INTERNAL_ENERGY_RELATIVE_FLOOR_J_M3", "2.3772267794722207e6"))
const CICE_ENTHALPY_RELATIVE_FLOOR =
    parse(Float64, get(ENV, "CICE_ENTHALPY_RELATIVE_FLOOR_J_M3", "2.958678993456372e8"))

struct NetCDFVariable
    dims :: Vector{String}
    shape :: Vector{Int}
end

struct NetCDFHeader
    dim_sizes :: Dict{String, Int}
    variables :: Dict{String, NetCDFVariable}
end

Base.@kwdef struct CICEHistoryRecord
    source_file :: String
    frequency :: String
    time_days :: Float64
    active_j :: Int
    active_i :: Int
    hi :: Float64
    hs :: Float64
    aice :: Float64
    Tsfc :: Float64
    congel :: Float64
    meltt :: Float64
    meltb :: Float64
    fsurf_ai :: Float64
    fcondtop_ai :: Float64
    fbot :: Float64
    siflcondbot :: Float64
    fhocn :: Float64
    fhocn_ai :: Float64
    fswabs :: Float64
    fswthru :: Float64
    fsens :: Float64
    flat :: Float64
    flwdn :: Float64
    flwup :: Float64
    Tair :: Float64
    sst :: Float64
    frzmlt :: Float64
    Tice :: Vector{Float64}
    Sice :: Vector{Float64}
    Tsnow :: Vector{Float64}
end

Base.@kwdef struct ReplayMetric
    metric :: String
    status :: String
    value :: Float64
    tolerance :: Float64
    evidence :: String
end

cice_case_id(history_dir) =
    CASE_ID == "" ? replace(basename(dirname(history_dir)), r"^bl99_|_(MU71|bubbly).*" => "") : CASE_ID

function cice_conductivity(::Type{FT}=Float64) where FT
    variant = lowercase(CONDUCTIVITY_VARIANT)
    if variant == "mu71"
        return MaykutUntersteinerConductivity(FT)
    elseif variant == "bubbly"
        return BubblyBrineConductivity(FT)
    else
        error("Unknown CICE_CONDUCT=$(CONDUCTIVITY_VARIANT). Expected MU71 or bubbly.")
    end
end

function ice_area_weight(record)
    return isfinite(record.aice) && record.aice > eps(Float64) ? record.aice : 1.0
end

function thermodynamic_ice_thickness(record)
    isfinite(record.hi) || return NaN
    return record.hi / ice_area_weight(record)
end

function per_ice_area_top_flux(record, flux)
    isfinite(flux) || return flux
    return flux / ice_area_weight(record)
end

function read_dashboard(path)
    isfile(path) || error("Dashboard not found: $path")
    lines = readlines(path)
    header = split(first(lines), ',')
    return [Dict(header .=> split(line, ','; limit=length(header)))
            for line in lines[2:end] if !isempty(line)]
end

function write_dashboard(path, rows)
    open(path, "w") do file
        println(file, "gate,status,evidence")
        for row in rows
            println(file, join((row["gate"], row["status"], row["evidence"]), ','))
        end
    end
end

function upsert_dashboard_row!(rows, new_row)
    index = findfirst(row -> row["gate"] == new_row["gate"], rows)
    if isnothing(index)
        comparison_index = findfirst(row -> row["gate"] == "cice_case_comparison", rows)
        insert!(rows, something(comparison_index, length(rows) + 1), new_row)
    else
        rows[index] = new_row
    end

    return rows
end

function print_dashboard(rows)
    println("Validation gates")
    println("================")
    for row in rows
        @printf("%-44s %-8s %s\n", row["gate"], row["status"], row["evidence"])
    end
end

function require_ncdump()
    success(pipeline(`which ncdump`; stdout=devnull, stderr=devnull)) ||
        error("`ncdump` is required to read CICE NetCDF history without adding package dependencies.")
end

function ncdump(args::AbstractString...)
    return read(Cmd(vcat(["ncdump"], collect(args))), String)
end

function parse_netcdf_header(path)
    output = ncdump("-h", path)
    dim_sizes = Dict{String, Int}()
    variables = Dict{String, NetCDFVariable}()

    for line in split(output, '\n')
        if (match_result = match(r"^\s*(\w+)\s*=\s*UNLIMITED\s*;\s*//\s*\((\d+)\s+currently\)", line)) !== nothing
            dim_sizes[match_result.captures[1]] = parse(Int, match_result.captures[2])
        elseif (match_result = match(r"^\s*(\w+)\s*=\s*(\d+)\s*;", line)) !== nothing
            dim_sizes[match_result.captures[1]] = parse(Int, match_result.captures[2])
        elseif (match_result = match(r"^\s*(?:byte|short|int|float|double)\s+(\w+)\(([^)]*)\)\s*;", line)) !== nothing
            name = match_result.captures[1]
            dims = strip.(split(match_result.captures[2], ','))
            variables[name] = NetCDFVariable(dims, [dim_sizes[dim] for dim in dims])
        end
    end

    return NetCDFHeader(dim_sizes, variables)
end

function parse_values(text)
    token_re = r"(?:_|[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?)"
    return [token.match == "_" ? NaN : parse(Float64, token.match)
            for token in eachmatch(token_re, text)]
end

function parse_ncdump_data(output)
    lines = split(output, '\n')
    data_start = findfirst(line -> strip(line) == "data:", lines)
    isnothing(data_start) && error("Could not find a data section in ncdump output.")

    blocks = Dict{String, Vector{Float64}}()
    current = nothing
    buffer = IOBuffer()

    function flush_current!()
        if current !== nothing
            blocks[current] = parse_values(String(take!(buffer)))
        end
    end

    for line in lines[(data_start + 1):end]
        stripped = strip(line)
        stripped == "}" && break

        if (match_result = match(r"^\s*([A-Za-z_]\w*)\s*=\s*(.*)$", line)) !== nothing
            flush_current!()
            current = match_result.captures[1]
            println(buffer, match_result.captures[2])
        elseif current !== nothing
            println(buffer, line)
        end
    end

    flush_current!()
    return blocks
end

function variable_blocks(path, variable_names)
    isempty(variable_names) && return Dict{String, Vector{Float64}}()
    return parse_ncdump_data(ncdump("-v", join(variable_names, ","), path))
end

function linear_position(shape, indices)
    length(shape) == length(indices) || error("shape and index lengths differ")
    position = 1
    stride = 1

    for dim in length(shape):-1:1
        index = indices[dim]
        1 <= index <= shape[dim] || error("Index $index outside dimension size $(shape[dim])")
        position += (index - 1) * stride
        stride *= shape[dim]
    end

    return position
end

function variable_value(header, blocks, name, index_by_dim)
    haskey(header.variables, name) || return NaN
    haskey(blocks, name) || return NaN

    variable = header.variables[name]
    indices = [get(index_by_dim, dim, 1) for dim in variable.dims]
    position = linear_position(variable.shape, indices)
    return blocks[name][position]
end

function active_cell(history_dir)
    grid_file = joinpath(history_dir, "iceh_grid.nc")
    isfile(grid_file) || error("CICE grid file not found: $grid_file")

    header = parse_netcdf_header(grid_file)
    blocks = variable_blocks(grid_file, ["tmask"])
    tmask = header.variables["tmask"]
    nj = header.dim_sizes["nj"]
    ni = header.dim_sizes["ni"]

    candidates = Tuple{Int, Int, Float64}[]
    center_j = (nj + 1) / 2
    center_i = (ni + 1) / 2

    for j in 1:tmask.shape[1], i in 1:tmask.shape[2]
        value = variable_value(header, blocks, "tmask", Dict("nj" => j, "ni" => i))
        if isfinite(value) && value == 1
            distance = (j - center_j)^2 + (i - center_i)^2
            push!(candidates, (j, i, distance))
        end
    end

    isempty(candidates) && error("No active T cell found in $grid_file")
    sort!(candidates; by=last)
    j, i, _ = first(candidates)
    return j, i
end

function available_names(header, names)
    return [name for name in names if haskey(header.variables, name)]
end

function scalar_field(header, blocks, name, j, i)
    haskey(header.variables, name) || return NaN
    return variable_value(header, blocks, name, Dict("time" => 1, "nj" => j, "ni" => i))
end

function profile_field(header, blocks, name, dim, j, i)
    haskey(header.variables, name) || return Float64[]
    length = header.dim_sizes[dim]

    return [variable_value(header, blocks, name,
                           Dict("time" => 1, "nc" => 1, dim => k, "nj" => j, "ni" => i))
            for k in 1:length]
end

function read_cice_record(path, j, i; frequency, suffix)
    header = parse_netcdf_header(path)
    candidate_names = ["time",
                       "hi$suffix", "hs$suffix", "aice$suffix", "Tsfc$suffix",
                       "congel$suffix", "meltt$suffix", "meltb$suffix",
                       "fsurf_ai$suffix", "fcondtop_ai$suffix",
                       "fbot$suffix", "siflcondbot$suffix",
                       "fhocn$suffix", "fhocn_ai$suffix",
                       "fswabs$suffix", "fswthru$suffix",
                       "fsens$suffix", "flat$suffix",
                       "flwdn$suffix", "flwup$suffix",
                       "Tair$suffix", "sst$suffix", "frzmlt$suffix",
                       "Tinz$suffix", "Sinz$suffix", "Tsnz$suffix"]
    names = available_names(header, candidate_names)
    blocks = variable_blocks(path, names)

    return CICEHistoryRecord(source_file=basename(path),
                             frequency=frequency,
                             time_days=variable_value(header, blocks, "time", Dict("time" => 1)),
                             active_j=j,
                             active_i=i,
                             hi=scalar_field(header, blocks, "hi$suffix", j, i),
                             hs=scalar_field(header, blocks, "hs$suffix", j, i),
                             aice=scalar_field(header, blocks, "aice$suffix", j, i),
                             Tsfc=scalar_field(header, blocks, "Tsfc$suffix", j, i),
                             congel=scalar_field(header, blocks, "congel$suffix", j, i),
                             meltt=scalar_field(header, blocks, "meltt$suffix", j, i),
                             meltb=scalar_field(header, blocks, "meltb$suffix", j, i),
                             fsurf_ai=scalar_field(header, blocks, "fsurf_ai$suffix", j, i),
                             fcondtop_ai=scalar_field(header, blocks, "fcondtop_ai$suffix", j, i),
                             fbot=scalar_field(header, blocks, "fbot$suffix", j, i),
                             siflcondbot=scalar_field(header, blocks, "siflcondbot$suffix", j, i),
                             fhocn=scalar_field(header, blocks, "fhocn$suffix", j, i),
                             fhocn_ai=scalar_field(header, blocks, "fhocn_ai$suffix", j, i),
                             fswabs=scalar_field(header, blocks, "fswabs$suffix", j, i),
                             fswthru=scalar_field(header, blocks, "fswthru$suffix", j, i),
                             fsens=scalar_field(header, blocks, "fsens$suffix", j, i),
                             flat=scalar_field(header, blocks, "flat$suffix", j, i),
                             flwdn=scalar_field(header, blocks, "flwdn$suffix", j, i),
                             flwup=scalar_field(header, blocks, "flwup$suffix", j, i),
                             Tair=scalar_field(header, blocks, "Tair$suffix", j, i),
                             sst=scalar_field(header, blocks, "sst$suffix", j, i),
                             frzmlt=scalar_field(header, blocks, "frzmlt$suffix", j, i),
                             Tice=profile_field(header, blocks, "Tinz$suffix", "nkice", j, i),
                             Sice=profile_field(header, blocks, "Sinz$suffix", "nkice", j, i),
                             Tsnow=profile_field(header, blocks, "Tsnz$suffix", "nksnow", j, i))
end

function discover_history_files(history_dir)
    isdir(history_dir) || return String[], String[]
    names = readdir(history_dir)
    ic_files = sort(joinpath.(Ref(history_dir),
                              filter(name -> startswith(name, "iceh_ic.") && endswith(name, ".nc"), names)))
    hourly_files = sort(joinpath.(Ref(history_dir),
                                  filter(name -> (startswith(name, "iceh_01h.") ||
                                                  startswith(name, "iceh_inst.")) &&
                                                 endswith(name, ".nc"), names)))
    return ic_files, hourly_files
end

function read_cice_history(history_dir)
    ic_files, hourly_files = discover_history_files(history_dir)
    isempty(ic_files) && error("No CICE initial-condition history file found in $history_dir")
    isempty(hourly_files) && error("No CICE hourly history files found in $history_dir")

    j, i = active_cell(history_dir)
    records = CICEHistoryRecord[]
    push!(records, read_cice_record(first(ic_files), j, i; frequency="initial", suffix=""))
    append!(records, [read_cice_record(path, j, i; frequency="hourly", suffix="_h") for path in hourly_files])
    sort!(records; by=record -> record.time_days)
    return records, j, i
end

function format_csv_value(value)
    return isnan(value) ? "" : @sprintf("%.12g", value)
end

function csv_field(value)
    text = string(value)
    if occursin(',', text) || occursin('"', text) || occursin('\n', text)
        return "\"" * replace(text, "\"" => "\"\"") * "\""
    else
        return text
    end
end

function write_cice_summary(path, records)
    mkpath(dirname(path))
    max_nice = maximum(length(record.Tice) for record in records)
    max_nsnow = maximum(length(record.Tsnow) for record in records)

    header = ["source_file", "frequency", "time_days", "active_j", "active_i",
              "hi_m", "hs_m", "aice", "thermodynamic_hi_m", "Tsfc_C", "congel_cm_per_day",
              "meltt_cm_per_day", "meltb_cm_per_day",
              "fsurf_ai_W_m2", "fcondtop_ai_W_m2",
              "fsurf_W_m2", "fcondtop_W_m2",
              "fbot_W_m2", "siflcondbot_W_m2",
              "fhocn_W_m2", "fhocn_ai_W_m2",
              "fswabs_W_m2", "fswthru_W_m2",
              "fsens_W_m2", "flat_W_m2",
              "flwdn_W_m2", "flwup_W_m2",
              "Tair_C", "sst_C", "frzmlt_W_m2"]
    append!(header, ["Tice_C_$k" for k in 1:max_nice])
    append!(header, ["Sice_ppt_$k" for k in 1:max_nice])
    append!(header, ["Tsnow_C_$k" for k in 1:max_nsnow])

    open(path, "w") do file
        println(file, join(header, ','))

        for record in records
            row = [record.source_file,
                   record.frequency,
                   format_csv_value(record.time_days),
                   string(record.active_j),
                   string(record.active_i),
                   format_csv_value(record.hi),
                   format_csv_value(record.hs),
                   format_csv_value(record.aice),
                   format_csv_value(thermodynamic_ice_thickness(record)),
                   format_csv_value(record.Tsfc),
                   format_csv_value(record.congel),
                   format_csv_value(record.meltt),
                   format_csv_value(record.meltb),
                   format_csv_value(record.fsurf_ai),
                   format_csv_value(record.fcondtop_ai),
                   format_csv_value(per_ice_area_top_flux(record, record.fsurf_ai)),
                   format_csv_value(per_ice_area_top_flux(record, record.fcondtop_ai)),
                   format_csv_value(record.fbot),
                   format_csv_value(record.siflcondbot),
                   format_csv_value(record.fhocn),
                   format_csv_value(record.fhocn_ai),
                   format_csv_value(record.fswabs),
                   format_csv_value(record.fswthru),
                   format_csv_value(record.fsens),
                   format_csv_value(record.flat),
                   format_csv_value(record.flwdn),
                   format_csv_value(record.flwup),
                   format_csv_value(record.Tair),
                   format_csv_value(record.sst),
                   format_csv_value(record.frzmlt)]

            append!(row, [k <= length(record.Tice) ? format_csv_value(record.Tice[k]) : "" for k in 1:max_nice])
            append!(row, [k <= length(record.Sice) ? format_csv_value(record.Sice[k]) : "" for k in 1:max_nice])
            append!(row, [k <= length(record.Tsnow) ? format_csv_value(record.Tsnow[k]) : "" for k in 1:max_nsnow])
            println(file, join(row, ','))
        end
    end
end

function finite_range(values)
    finite_values = filter(isfinite, values)
    isempty(finite_values) && return (NaN, NaN)
    return extrema(finite_values)
end

function finite_maximum(values)
    finite_values = filter(isfinite, values)
    isempty(finite_values) && return NaN
    return maximum(finite_values)
end

function finite_mean(values)
    finite_values = filter(isfinite, values)
    isempty(finite_values) && return NaN
    return sum(finite_values) / length(finite_values)
end

function cice_thickness_rate_integrals(records)
    integrated_basal_growth = 0.0
    integrated_top_melt = 0.0
    integrated_bottom_melt = 0.0

    for (start_record, target_record) in zip(records[1:end-1], records[2:end])
        Δt_days = target_record.time_days - start_record.time_days
        isfinite(Δt_days) && Δt_days > 0 || continue

        if isfinite(target_record.congel)
            integrated_basal_growth += target_record.congel * Δt_days / 100
        end

        if isfinite(target_record.meltt)
            integrated_top_melt += target_record.meltt * Δt_days / 100
        end

        if isfinite(target_record.meltb)
            integrated_bottom_melt += target_record.meltb * Δt_days / 100
        end
    end

    rate_integrated_delta =
        integrated_basal_growth - integrated_top_melt - integrated_bottom_melt
    observed_delta = last(records).hi - first(records).hi
    residual = observed_delta - rate_integrated_delta

    return (; integrated_basal_growth,
            integrated_top_melt,
            integrated_bottom_melt,
            rate_integrated_delta,
            observed_delta,
            residual)
end

cice_thickness_budget_status(records) =
    abs(cice_thickness_rate_integrals(records).residual) <= 1e-3 ? "pass" : "fail"

function thickness_acceptance_ratio(predicted, cice)
    absolute_error = abs(predicted - cice)

    if abs(cice) < 0.1
        return absolute_error / 1e-3
    else
        return absolute_error / (0.01 * abs(cice))
    end
end

strict_bl99_validation_status(case_id, forced_replay_metrics) =
    metric_status(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_acceptance_ratio") == "pass" &&
    metric_status(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_temperature_acceptance_ratio") == "pass" &&
    metric_status(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_cice_enthalpy_error") == "pass" &&
    metric_status(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_column_energy_error") == "pass" ?
    "pass" : "incomplete"

function print_cice_summary(records, history_dir, summary_path)
    first_record = first(records)
    final_record = last(records)
    initial_tmin, initial_tmax = finite_range(first_record.Tice)
    final_tmin, final_tmax = finite_range(final_record.Tice)
    max_growth = finite_maximum([record.congel for record in records])
    max_top_melt = finite_maximum([record.meltt for record in records])
    max_bottom_melt = finite_maximum([record.meltb for record in records])
    thickness_rates = cice_thickness_rate_integrals(records)
    first_hourly = records[findfirst(record -> record.frequency == "hourly", records)]

    println()
    println("CICE history parser")
    println("===================")
    println("history_dir: $history_dir")
    println("summary_csv: $summary_path")
    @printf("records: %d (%d hourly)\n", length(records), count(record -> record.frequency == "hourly", records))
    @printf("active T cell: j=%d i=%d\n", first_record.active_j, first_record.active_i)
    @printf("initial: t=%.8g d hi=%.8g m aice=%.8g thermodynamic_hi=%.8g m hs=%.8g m Tsfc=%.8g C Tice=[%.8g, %.8g] C\n",
            first_record.time_days, first_record.hi, first_record.aice,
            thermodynamic_ice_thickness(first_record), first_record.hs, first_record.Tsfc,
            initial_tmin, initial_tmax)
    @printf("final:   t=%.8g d hi=%.8g m aice=%.8g thermodynamic_hi=%.8g m hs=%.8g m Tsfc=%.8g C Tice=[%.8g, %.8g] C\n",
            final_record.time_days, final_record.hi, final_record.aice,
            thermodynamic_ice_thickness(final_record), final_record.hs, final_record.Tsfc,
            final_tmin, final_tmax)
    @printf("max rates: congel=%.8g cm/day meltt=%.8g cm/day meltb=%.8g cm/day\n",
            max_growth, max_top_melt, max_bottom_melt)
    @printf("rate-integrated thickness: congel=%.8g m meltt=%.8g m meltb=%.8g m net=%.8g m residual=%.8g m\n",
            thickness_rates.integrated_basal_growth,
            thickness_rates.integrated_top_melt,
            thickness_rates.integrated_bottom_melt,
            thickness_rates.rate_integrated_delta,
            thickness_rates.residual)
    @printf("first-hour fluxes: fsurf_ai=%.8g W/m^2 fcondtop_ai=%.8g W/m^2 aice=%.8g fsurf=%.8g W/m^2 fcondtop=%.8g W/m^2 fhocn=%.8g W/m^2\n",
            first_hourly.fsurf_ai,
            first_hourly.fcondtop_ai,
            first_hourly.aice,
            per_ice_area_top_flux(first_hourly, first_hourly.fsurf_ai),
            per_ice_area_top_flux(first_hourly, first_hourly.fcondtop_ai),
            first_hourly.fhocn)
end

function cice_profile_function(values_top_to_bottom; thickness=1)
    values = copy(values_top_to_bottom)
    n = length(values)

    return z -> begin
        zeta = clamp(z / thickness, 0, 1)
        depth_from_surface = 1 - zeta
        top_index = clamp(floor(Int, depth_from_surface * n) + 1, 1, n)
        values[top_index]
    end
end

column_values(field) = vec(Array(interior(field, 1, 1, :)))

function top_to_bottom_profile(field)
    return reverse(column_values(field))
end

function max_salinity_profile_error(records)
    profile = FixedDrainedIceSalinityProfile(Float64)
    max_error = 0.0

    for record in records
        n = length(record.Sice)
        for k in 1:n
            normalized_depth = (k - 0.5) / n
            expected = salinity_at_normalized_depth(profile, normalized_depth)
            max_error = max(max_error, abs(record.Sice[k] - expected))
        end
    end

    return max_error
end

function cice_internal_energy_profile(relation, record)
    return [internal_energy(relation, record.Tice[k], record.Sice[k])
            for k in eachindex(record.Tice)]
end

function predicted_internal_energy_profile(relation, predicted_temperature, target_record)
    return [internal_energy(relation, predicted_temperature[k], target_record.Sice[k])
            for k in eachindex(predicted_temperature)]
end

function internal_energy_relative_floor(relation, records)
    if isfinite(INTERNAL_ENERGY_RELATIVE_FLOOR) && INTERNAL_ENERGY_RELATIVE_FLOOR > 0
        return INTERNAL_ENERGY_RELATIVE_FLOOR
    end

    energies = reduce(vcat, [abs.(cice_internal_energy_profile(relation, record))
                             for record in records])
    nonzero_energies = filter(energy -> isfinite(energy) && energy > 0, energies)
    return isempty(nonzero_energies) ? 1.0 : minimum(nonzero_energies)
end

column_integrated_energy(energy_top_to_bottom, thickness) =
    sum(energy_top_to_bottom) * thickness / length(energy_top_to_bottom)

function cice_enthalpy_offset(relation)
    phase_transitions = relation.phase_transitions
    return phase_transitions.density * phase_transitions.reference_latent_heat
end

function replay_profile_errors(predicted_temperature,
                               target_record,
                               relation,
                               energy_floor;
                               enthalpy_floor = CICE_ENTHALPY_RELATIVE_FLOOR)
    cice_energy = cice_internal_energy_profile(relation, target_record)
    predicted_energy =
        predicted_internal_energy_profile(relation, predicted_temperature, target_record)

    absolute_temperature_errors = abs.(predicted_temperature .- target_record.Tice)
    relative_temperature_errors = absolute_temperature_errors ./ max.(abs.(target_record.Tice), 1)
    absolute_temperature_acceptance_ratios = absolute_temperature_errors ./ 0.02
    relative_temperature_acceptance_ratios = relative_temperature_errors ./ 0.01
    temperature_acceptance_ratios =
        [abs(target_record.Tice[k]) < 2 ?
         min(relative_temperature_acceptance_ratios[k], absolute_temperature_acceptance_ratios[k]) :
         relative_temperature_acceptance_ratios[k]
         for k in eachindex(target_record.Tice)]

    absolute_energy_errors = abs.(predicted_energy .- cice_energy)
    relative_energy_errors = absolute_energy_errors ./ max.(abs.(cice_energy), energy_floor)
    enthalpy_offset = cice_enthalpy_offset(relation)
    cice_enthalpy = cice_energy .- enthalpy_offset
    predicted_enthalpy = predicted_energy .- enthalpy_offset
    absolute_enthalpy_errors = abs.(predicted_enthalpy .- cice_enthalpy)
    relative_enthalpy_errors = absolute_enthalpy_errors ./ max.(abs.(cice_enthalpy), enthalpy_floor)

    cice_column_energy = column_integrated_energy(cice_energy, target_record.hi)
    predicted_column_energy = column_integrated_energy(predicted_energy, target_record.hi)
    column_energy_floor = energy_floor * target_record.hi
    relative_column_energy_error =
        abs(predicted_column_energy - cice_column_energy) /
        max(abs(cice_column_energy), column_energy_floor)

    return (; absolute_temperature_errors,
            relative_temperature_errors,
            temperature_acceptance_ratios,
            absolute_energy_errors,
            relative_energy_errors,
            absolute_enthalpy_errors,
            relative_enthalpy_errors,
            relative_column_energy_error)
end

function run_climaseaice_state_replay(records)
    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    first_record = first(records)
    n = length(first_record.Tice)
    grid = RectilinearGrid(size = n,
                           z = (0, 1),
                           topology = (Flat, Flat, Bounded))

    energy_transport =
        ConductiveTemperatureTransport(conductivity = cice_conductivity(Float64))

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = 0.0,
        energy_transport)

    max_ingest_temperature_error = 0.0
    max_ingest_salinity_error = 0.0
    max_relation_roundtrip_error = 0.0
    max_reconstructed_energy_error = 0.0

    for record in records
        set!(thermodynamics;
             bulk_salinity = cice_profile_function(record.Sice),
             temperature = cice_profile_function(record.Tice))

        compute_column_transport_coefficients!(thermodynamics)

        climat_t = column_values(thermodynamics.fields.temperature)
        climat_s = column_values(thermodynamics.fields.bulk_salinity)
        climat_e = column_values(thermodynamics.fields.internal_energy)
        cice_t_bottom_to_top = reverse(record.Tice)
        cice_s_bottom_to_top = reverse(record.Sice)

        max_ingest_temperature_error =
            max(max_ingest_temperature_error,
                maximum(abs.(climat_t .- cice_t_bottom_to_top)))

        max_ingest_salinity_error =
            max(max_ingest_salinity_error,
                maximum(abs.(climat_s .- cice_s_bottom_to_top)))

        for k in eachindex(climat_t)
            E = internal_energy(relation, climat_t[k], climat_s[k])
            cice_energy = internal_energy(relation, cice_t_bottom_to_top[k], cice_s_bottom_to_top[k])
            recovered_temperature = temperature(relation, E, climat_s[k])
            max_relation_roundtrip_error =
                max(max_relation_roundtrip_error,
                    abs(recovered_temperature - climat_t[k]),
                    abs(temperature(relation, climat_e[k], climat_s[k]) - climat_t[k]))
            max_reconstructed_energy_error =
                max(max_reconstructed_energy_error, abs(E - cice_energy))
        end
    end

    metrics = ReplayMetric[
        ReplayMetric(metric = "cice_salinity_profile_max_abs_error_ppt",
                     status = max_salinity_profile_error(records) <= 1e-5 ? "pass" : "fail",
                     value = max_salinity_profile_error(records),
                     tolerance = 1e-5,
                     evidence = "CICE Sice profile vs FixedDrainedIceSalinityProfile layer midpoints"),
        ReplayMetric(metric = "climaseaice_column_ingest_temperature_max_abs_error_C",
                     status = max_ingest_temperature_error <= 1e-12 ? "pass" : "fail",
                     value = max_ingest_temperature_error,
                     tolerance = 1e-12,
                     evidence = "CICE Tinz mapped top-to-bottom into Oceananigans bottom-to-top column"),
        ReplayMetric(metric = "climaseaice_column_ingest_salinity_max_abs_error_ppt",
                     status = max_ingest_salinity_error <= 1e-12 ? "pass" : "fail",
                     value = max_ingest_salinity_error,
                     tolerance = 1e-12,
                     evidence = "CICE Sinz mapped top-to-bottom into Oceananigans bottom-to-top column"),
        ReplayMetric(metric = "fixed_relation_temperature_roundtrip_max_abs_error_C",
                     status = max_relation_roundtrip_error <= 1e-11 ? "pass" : "fail",
                     value = max_relation_roundtrip_error,
                     tolerance = 1e-11,
                     evidence = "FixedSalinityBrinePocketEnergyRelation energy-temperature round trip"),
        ReplayMetric(metric = "reconstructed_internal_energy_max_abs_error_J_m3",
                     status = max_reconstructed_energy_error <= 1e-7 ? "pass" : "fail",
                     value = max_reconstructed_energy_error,
                     tolerance = 1e-7,
                     evidence = "CICE Tinz/Sinz reconstructed through FixedSalinityBrinePocketEnergyRelation"),
    ]

    return metrics
end

function required_fluxes_available(records)
    hourly_records = filter(record -> record.frequency == "hourly", records)
    return all(record -> isfinite(record.fcondtop_ai), hourly_records)
end

function cice_bottom_conductive_flux(start_record,
                                     target_record,
                                     conductivity = cice_conductivity(Float64))
    n = length(start_record.Tice)
    hilyr = start_record.hi / n
    bottom_conductivity = ice_thermal_conductivity(conductivity,
                                                   last(start_record.Tice),
                                                   last(start_record.Sice))
    bottom_conductance = 2 * bottom_conductivity / hilyr
    Tbot = isfinite(target_record.sst) ? target_record.sst : start_record.sst
    return bottom_conductance * (last(target_record.Tice) - Tbot)
end

function icepack_temperature_matrix_replay_step(start_record,
                                                target_record,
                                                relation,
                                                energy_transport;
                                                start_temperature = start_record.Tice,
                                                remap_to_target_thickness = false)
    Δt = (target_record.time_days - start_record.time_days) * 86400
    thickness = start_record.hi
    Tbot = isfinite(target_record.sst) ? target_record.sst : start_record.sst
    grid = RectilinearGrid(size = length(start_temperature),
                           z = (0, thickness),
                           topology = (Flat, Flat, Bounded))
    boundary_conditions =
        ColumnBoundaryConditions(top = PrescribedEnergyFlux(flux = per_ice_area_top_flux(target_record,
                                                                                         target_record.fcondtop_ai)),
                                 bottom = PrescribedTemperature(Tbot))
    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = cice_profile_function(start_record.Sice; thickness),
        energy_transport,
        boundary_conditions)

    set!(thermodynamics;
         bulk_salinity = cice_profile_function(start_record.Sice; thickness),
         temperature = cice_profile_function(start_temperature; thickness))

    icepack_temperature_matrix_step!(thermodynamics, Δt; tolerance = 1e-12)

    if remap_to_target_thickness && !isapprox(target_record.hi, start_record.hi; atol = 1e-12, rtol = 0)
        return remap_temperature_to_target_thickness(thermodynamics,
                                                    start_record,
                                                    target_record,
                                                    relation)
    else
        return top_to_bottom_profile(thermodynamics.fields.temperature)
    end
end

layer_faces(thickness, n; offset = 0) =
    collect(range(offset, offset + thickness; length = n + 1))

function basal_growth_fill_energy(relation, target_record)
    # Icepack BL99 basal congelation adds ice at Tbot with bottom-layer salinity,
    # then repartitions to equal layers through adjust_enthalpy.
    Tbot = isfinite(target_record.sst) ? target_record.sst : last(target_record.Tice)
    return internal_energy(relation, Tbot, last(target_record.Sice))
end

function remap_temperature_to_target_thickness(thermodynamics,
                                               start_record,
                                               target_record,
                                               relation)
    source_energy = column_values(thermodynamics.fields.internal_energy)
    n = length(source_energy)
    source_thickness = start_record.hi
    target_thickness = target_record.hi
    thickness_change = target_thickness - source_thickness
    target_faces = layer_faces(target_thickness, n)

    if thickness_change > 0
        source_faces = layer_faces(source_thickness, n; offset = thickness_change)
        fill_energy = basal_growth_fill_energy(relation, target_record)
    else
        source_faces = layer_faces(source_thickness, n)
        fill_energy = zero(eltype(source_energy))
    end

    column_energy_thickness_remap!(thermodynamics,
                                   source_faces,
                                   target_faces;
                                   fill_energy,
                                   bulk_salinity = cice_profile_function(target_record.Sice;
                                                                         thickness = source_thickness))

    return top_to_bottom_profile(thermodynamics.fields.temperature)
end

function set_column_metric!(grid, previous_metric, current_metric; surface = 0)
    fill!(grid.z.ηⁿ, surface)
    fill!(grid.z.σᶜᶜ⁻, previous_metric)
    fill!(grid.z.σᶜᶜⁿ, current_metric)
    fill!(grid.z.σᶠᶜⁿ, current_metric)
    fill!(grid.z.σᶜᶠⁿ, current_metric)
    fill!(grid.z.σᶠᶠⁿ, current_metric)
    return nothing
end

function moving_metric_forced_step(start_record,
                                   start_temperature,
                                   target_record,
                                   relation,
                                   energy_transport)
    Δt = (target_record.time_days - start_record.time_days) * 86400
    reference_thickness = start_record.hi
    target_metric = target_record.hi / reference_thickness
    Tbot = isfinite(target_record.sst) ? target_record.sst : start_record.sst
    grid = RectilinearGrid(size = length(start_temperature),
                           z = MutableVerticalDiscretization((0, reference_thickness)),
                           topology = (Flat, Flat, Bounded))

    boundary_conditions =
        ColumnBoundaryConditions(top = PrescribedEnergyFlux(flux = per_ice_area_top_flux(target_record,
                                                                                         target_record.fcondtop_ai)),
                                 bottom = PrescribedTemperature(Tbot))

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = cice_profile_function(start_record.Sice; thickness = reference_thickness),
        energy_transport,
        boundary_conditions)

    set_column_metric!(grid, 1.0, 1.0)
    set!(thermodynamics;
         bulk_salinity = cice_profile_function(start_record.Sice; thickness = reference_thickness),
         temperature = cice_profile_function(start_temperature; thickness = reference_thickness))

    set_column_metric!(grid, 1.0, target_metric)
    column_energy_time_step!(thermodynamics, Δt)

    return top_to_bottom_profile(thermodynamics.fields.temperature)
end

function predicted_bottom_conductive_flux(temperature_top_to_bottom,
                                          salinity_top_to_bottom,
                                          thickness,
                                          target_record,
                                          conductivity = cice_conductivity(Float64))
    n = length(temperature_top_to_bottom)
    hilyr = thickness / n
    bottom_conductivity = ice_thermal_conductivity(conductivity,
                                                   last(temperature_top_to_bottom),
                                                   last(salinity_top_to_bottom))
    bottom_conductance = 2 * bottom_conductivity / hilyr
    Tbot = isfinite(target_record.sst) ? target_record.sst : last(temperature_top_to_bottom)
    return bottom_conductance * (last(temperature_top_to_bottom) - Tbot)
end

function cice_negative_enthalpy_magnitude(relation, temperature, salinity)
    shifted_energy = internal_energy(relation, temperature, salinity)
    negative_enthalpy = shifted_energy - cice_enthalpy_offset(relation)
    return max(eps(abs(negative_enthalpy)), -negative_enthalpy)
end

function prognostic_thickness_forced_step(start_record,
                                          start_temperature,
                                          start_thickness,
                                          target_record,
                                          relation,
                                          energy_transport;
                                          fixed_thickness = false)
    Δt = (target_record.time_days - start_record.time_days) * 86400
    requested_surface_flux = isfinite(target_record.fsurf_ai) ?
                             per_ice_area_top_flux(target_record, target_record.fsurf_ai) :
                             per_ice_area_top_flux(target_record, target_record.fcondtop_ai)
    applied_surface_flux = isfinite(target_record.fcondtop_ai) ?
                           per_ice_area_top_flux(target_record, target_record.fcondtop_ai) :
                           requested_surface_flux
    Tbot = isfinite(target_record.sst) ? target_record.sst : start_record.sst
    grid = RectilinearGrid(size = length(start_temperature),
                           z = (0, start_thickness),
                           topology = (Flat, Flat, Bounded))

    boundary_conditions =
        ColumnBoundaryConditions(top = PrescribedEnergyFlux(flux = applied_surface_flux),
                                 bottom = PrescribedTemperature(Tbot))

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = cice_profile_function(start_record.Sice; thickness = start_thickness),
        energy_transport,
        boundary_conditions)

    set!(thermodynamics;
         bulk_salinity = cice_profile_function(start_record.Sice; thickness = start_thickness),
         temperature = cice_profile_function(start_temperature; thickness = start_thickness))

    surface_residual_flux = applied_surface_flux - requested_surface_flux
    ocean_heat_flux = isfinite(target_record.fhocn) ? target_record.fhocn : 0.0

    column_energy_time_step!(thermodynamics, Δt)
    stepped_temperature = top_to_bottom_profile(thermodynamics.fields.temperature)
    stepped_salinity = top_to_bottom_profile(thermodynamics.fields.bulk_salinity)
    bottom_conductive_flux =
        predicted_bottom_conductive_flux(stepped_temperature,
                                         stepped_salinity,
                                         start_thickness,
                                         target_record)
    basal_residual_flux = -bottom_conductive_flux + ocean_heat_flux

    surface_enthalpy_magnitude =
        cice_negative_enthalpy_magnitude(relation, first(stepped_temperature), first(stepped_salinity))
    basal_temperature = if basal_residual_flux >= 0
        isfinite(target_record.sst) ? target_record.sst : last(start_temperature)
    else
        last(stepped_temperature)
    end
    basal_enthalpy_magnitude =
        cice_negative_enthalpy_magnitude(relation, basal_temperature, last(stepped_salinity))

    surface_thickness_change = fixed_thickness ? 0.0 :
        Δt * surface_residual_flux / surface_enthalpy_magnitude

    basal_thickness_change = fixed_thickness ? 0.0 :
        Δt * basal_residual_flux / basal_enthalpy_magnitude

    target_thickness = max(start_thickness + surface_thickness_change + basal_thickness_change, 1e-6)

    if !isapprox(target_thickness, start_thickness; atol = 1e-12, rtol = 0)
        n = length(start_temperature)
        source_faces = layer_faces(start_thickness, n; offset = basal_thickness_change)
        target_faces = layer_faces(target_thickness, n)
        fill_energy = basal_thickness_change > 0 ?
                      basal_growth_fill_energy(relation, target_record) :
                      zero(eltype(grid))

        column_energy_thickness_remap!(thermodynamics,
                                       source_faces,
                                       target_faces;
                                       fill_energy,
                                       bulk_salinity = cice_profile_function(target_record.Sice;
                                                                             thickness = start_thickness))
    end

    return (temperature = top_to_bottom_profile(thermodynamics.fields.temperature),
            thickness = target_thickness,
            surface_thickness_change,
            basal_thickness_change,
            surface_residual_flux,
            basal_residual_flux)
end

function icepack_matrix_prognostic_thickness_forced_step(start_record,
                                                         start_temperature,
                                                         start_thickness,
                                                         target_record,
                                                         relation,
                                                         energy_transport;
                                                         fixed_thickness = false)
    Δt = (target_record.time_days - start_record.time_days) * 86400
    requested_surface_flux = isfinite(target_record.fsurf_ai) ?
                             per_ice_area_top_flux(target_record, target_record.fsurf_ai) :
                             per_ice_area_top_flux(target_record, target_record.fcondtop_ai)
    applied_surface_flux = isfinite(target_record.fcondtop_ai) ?
                           per_ice_area_top_flux(target_record, target_record.fcondtop_ai) :
                           requested_surface_flux
    Tbot = isfinite(target_record.sst) ? target_record.sst : start_record.sst
    grid = RectilinearGrid(size = length(start_temperature),
                           z = (0, start_thickness),
                           topology = (Flat, Flat, Bounded))

    boundary_conditions =
        ColumnBoundaryConditions(top = PrescribedEnergyFlux(flux = applied_surface_flux),
                                 bottom = PrescribedTemperature(Tbot))

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = cice_profile_function(start_record.Sice; thickness = start_thickness),
        energy_transport,
        boundary_conditions)

    set!(thermodynamics;
         bulk_salinity = cice_profile_function(start_record.Sice; thickness = start_thickness),
         temperature = cice_profile_function(start_temperature; thickness = start_thickness))

    icepack_temperature_matrix_step!(thermodynamics, Δt; tolerance = 1e-12)

    stepped_temperature = top_to_bottom_profile(thermodynamics.fields.temperature)
    stepped_salinity = top_to_bottom_profile(thermodynamics.fields.bulk_salinity)
    surface_residual_flux = applied_surface_flux - requested_surface_flux
    bottom_conductive_flux =
        predicted_bottom_conductive_flux(stepped_temperature,
                                         stepped_salinity,
                                         start_thickness,
                                         target_record)
    ocean_heat_flux = isfinite(target_record.fhocn) ? target_record.fhocn : 0.0
    basal_residual_flux = -bottom_conductive_flux + ocean_heat_flux

    surface_enthalpy_magnitude =
        cice_negative_enthalpy_magnitude(relation, first(stepped_temperature), first(stepped_salinity))
    basal_temperature = basal_residual_flux >= 0 ? Tbot : last(stepped_temperature)
    basal_enthalpy_magnitude =
        cice_negative_enthalpy_magnitude(relation, basal_temperature, last(stepped_salinity))

    surface_thickness_change = fixed_thickness ? 0.0 :
        Δt * surface_residual_flux / surface_enthalpy_magnitude

    basal_thickness_change = fixed_thickness ? 0.0 :
        Δt * basal_residual_flux / basal_enthalpy_magnitude

    target_thickness = max(start_thickness + surface_thickness_change + basal_thickness_change, 1e-6)

    if !isapprox(target_thickness, start_thickness; atol = 1e-12, rtol = 0)
        n = length(start_temperature)
        source_faces = layer_faces(start_thickness, n; offset = basal_thickness_change)
        target_faces = layer_faces(target_thickness, n)
        fill_energy = basal_thickness_change > 0 ?
                      basal_growth_fill_energy(relation, target_record) :
                      zero(eltype(grid))

        column_energy_thickness_remap!(thermodynamics,
                                       source_faces,
                                       target_faces;
                                       fill_energy,
                                       bulk_salinity = cice_profile_function(target_record.Sice;
                                                                             thickness = start_thickness))
    end

    return (temperature = top_to_bottom_profile(thermodynamics.fields.temperature),
            thickness = target_thickness,
            surface_thickness_change,
            basal_thickness_change,
            surface_residual_flux,
            basal_residual_flux)
end

function fixed_grid_forced_step(start_record,
                                start_temperature,
                                target_record,
                                relation,
                                energy_transport;
                                remap_to_target_thickness = false)
    Δt = (target_record.time_days - start_record.time_days) * 86400
    thickness = start_record.hi
    Tbot = isfinite(target_record.sst) ? target_record.sst : start_record.sst
    grid = RectilinearGrid(size = length(start_temperature),
                           z = (0, thickness),
                           topology = (Flat, Flat, Bounded))

    boundary_conditions =
        ColumnBoundaryConditions(top = PrescribedEnergyFlux(flux = per_ice_area_top_flux(target_record,
                                                                                         target_record.fcondtop_ai)),
                                 bottom = PrescribedTemperature(Tbot))

    thermodynamics = prescribed_salinity_enthalpy_thermodynamics(grid;
        relation,
        salinity_profile = cice_profile_function(start_record.Sice; thickness),
        energy_transport,
        boundary_conditions)

    set!(thermodynamics;
         bulk_salinity = cice_profile_function(start_record.Sice; thickness),
         temperature = cice_profile_function(start_temperature; thickness))

    column_energy_time_step!(thermodynamics, Δt)

    if remap_to_target_thickness && !isapprox(target_record.hi, start_record.hi; atol = 1e-12, rtol = 0)
        return remap_temperature_to_target_thickness(thermodynamics,
                                                    start_record,
                                                    target_record,
                                                    relation)
    else
        return top_to_bottom_profile(thermodynamics.fields.temperature)
    end
end

function run_fixed_grid_forced_replay(records)
    hourly_records = filter(record -> record.frequency == "hourly", records)
    initial_record = first(records)

    if isempty(hourly_records) || !required_fluxes_available(records)
        return ReplayMetric[
            ReplayMetric(metric = "fixed_grid_forced_replay_available",
                         status = "not_run",
                         value = 0,
                         tolerance = 0,
                         evidence = "Missing hourly fcondtop_ai forcing fields"),
        ]
    end

    relation = FixedSalinityBrinePocketEnergyRelation(Float64)
    energy_floor = internal_energy_relative_floor(relation, records)
    energy_transport =
        ConductiveTemperatureTransport(conductivity = cice_conductivity(Float64))
    fixed_thickness_case = CASE_ID == "cold_conductive_relaxation"

    current_record = initial_record
    current_temperature = copy(initial_record.Tice)
    moving_record = initial_record
    moving_temperature = copy(initial_record.Tice)
    icepack_matrix_moving_record = initial_record
    icepack_matrix_moving_temperature = copy(initial_record.Tice)
    moving_metric_record = initial_record
    moving_metric_temperature = copy(initial_record.Tice)
    prognostic_record = initial_record
    prognostic_temperature = copy(initial_record.Tice)
    prognostic_thickness = initial_record.hi
    prognostic_top_melt = 0.0
    prognostic_basal_growth = 0.0
    prognostic_bottom_melt = 0.0
    prognostic_max_temperature_acceptance_ratio = 0.0
    prognostic_max_relative_cice_enthalpy_error = 0.0
    prognostic_max_relative_column_energy_error = 0.0
    icepack_matrix_prognostic_record = initial_record
    icepack_matrix_prognostic_temperature = copy(initial_record.Tice)
    icepack_matrix_prognostic_thickness = initial_record.hi
    icepack_matrix_prognostic_top_melt = 0.0
    icepack_matrix_prognostic_basal_growth = 0.0
    icepack_matrix_prognostic_bottom_melt = 0.0
    icepack_matrix_prognostic_max_temperature_acceptance_ratio = 0.0
    icepack_matrix_prognostic_max_relative_cice_enthalpy_error = 0.0
    icepack_matrix_prognostic_max_relative_column_energy_error = 0.0
    max_abs_temperature_error = 0.0
    max_relative_temperature_error = 0.0
    max_temperature_acceptance_ratio = 0.0
    max_relative_internal_energy_error = 0.0
    max_relative_cice_enthalpy_error = 0.0
    max_relative_column_energy_error = 0.0
    final_abs_temperature_error = 0.0
    final_relative_temperature_error = 0.0
    moving_max_abs_temperature_error = 0.0
    moving_max_relative_temperature_error = 0.0
    moving_max_temperature_acceptance_ratio = 0.0
    moving_max_relative_internal_energy_error = 0.0
    moving_max_relative_cice_enthalpy_error = 0.0
    moving_max_relative_column_energy_error = 0.0
    moving_final_abs_temperature_error = 0.0
    moving_final_relative_temperature_error = 0.0
    moving_metric_max_abs_temperature_error = 0.0
    moving_metric_max_relative_temperature_error = 0.0
    moving_metric_max_temperature_acceptance_ratio = 0.0
    moving_metric_max_relative_internal_energy_error = 0.0
    moving_metric_max_relative_cice_enthalpy_error = 0.0
    moving_metric_max_relative_column_energy_error = 0.0
    moving_metric_final_abs_temperature_error = 0.0
    moving_metric_final_relative_temperature_error = 0.0
    icepack_matrix_moving_max_abs_temperature_error = 0.0
    icepack_matrix_moving_max_relative_temperature_error = 0.0
    icepack_matrix_moving_max_relative_internal_energy_error = 0.0
    icepack_matrix_moving_max_relative_column_energy_error = 0.0
    icepack_matrix_moving_max_temperature_acceptance_ratio = 0.0
    icepack_matrix_moving_final_abs_temperature_error = 0.0
    icepack_matrix_moving_final_relative_temperature_error = 0.0
    icepack_matrix_moving_final_temperature_acceptance_ratio = 0.0
    icepack_matrix_moving_max_relative_cice_enthalpy_error = 0.0
    icepack_matrix_moving_max_error_time_days = NaN
    icepack_matrix_moving_max_error_layer_from_top = NaN
    icepack_matrix_moving_max_error_predicted_temperature = NaN
    icepack_matrix_moving_max_error_cice_temperature = NaN
    icepack_matrix_moving_final_top_temperature_error = NaN
    icepack_matrix_moving_final_bottom_temperature_error = NaN
    reset_max_abs_temperature_error = 0.0
    reset_max_relative_temperature_error = 0.0
    reset_max_temperature_acceptance_ratio = 0.0
    reset_max_relative_internal_energy_error = 0.0
    reset_max_relative_cice_enthalpy_error = 0.0
    reset_max_relative_column_energy_error = 0.0
    reset_first_abs_temperature_error = 0.0
    reset_first_relative_temperature_error = 0.0
    icepack_matrix_max_abs_temperature_error = 0.0
    icepack_matrix_max_relative_temperature_error = 0.0
    icepack_matrix_max_temperature_acceptance_ratio = 0.0
    icepack_matrix_max_relative_internal_energy_error = 0.0
    icepack_matrix_max_relative_cice_enthalpy_error = 0.0
    icepack_matrix_max_relative_column_energy_error = 0.0
    icepack_matrix_remapped_reset_max_abs_temperature_error = 0.0
    icepack_matrix_remapped_reset_max_relative_temperature_error = 0.0
    icepack_matrix_remapped_reset_max_temperature_acceptance_ratio = 0.0
    icepack_matrix_remapped_reset_max_relative_internal_energy_error = 0.0
    icepack_matrix_remapped_reset_max_relative_cice_enthalpy_error = 0.0
    icepack_matrix_remapped_reset_max_relative_column_energy_error = 0.0
    icepack_matrix_remapped_reset_intervals = 0
    thickness_remapped_reset_max_abs_temperature_error = 0.0
    thickness_remapped_reset_max_relative_temperature_error = 0.0
    thickness_remapped_reset_max_temperature_acceptance_ratio = 0.0
    thickness_remapped_reset_max_relative_internal_energy_error = 0.0
    thickness_remapped_reset_max_relative_cice_enthalpy_error = 0.0
    thickness_remapped_reset_max_relative_column_energy_error = 0.0
    thickness_remapped_reset_intervals = 0
    remapped_reset_max_abs_temperature_error = 0.0
    remapped_reset_max_relative_temperature_error = 0.0
    remapped_reset_max_temperature_acceptance_ratio = 0.0
    remapped_reset_max_relative_internal_energy_error = 0.0
    remapped_reset_max_relative_cice_enthalpy_error = 0.0
    remapped_reset_max_relative_column_energy_error = 0.0
    remapped_reset_intervals = 0
    first_bottom_conductive_flux = cice_bottom_conductive_flux(initial_record, first(hourly_records))

    for target_record in hourly_records
        predicted_temperature =
            fixed_grid_forced_step(current_record, current_temperature,
                                   target_record, relation, energy_transport)

        errors = replay_profile_errors(predicted_temperature, target_record, relation, energy_floor)
        final_abs_temperature_error = maximum(errors.absolute_temperature_errors)
        final_relative_temperature_error = maximum(errors.relative_temperature_errors)
        max_abs_temperature_error = max(max_abs_temperature_error, final_abs_temperature_error)
        max_relative_temperature_error = max(max_relative_temperature_error, final_relative_temperature_error)
        max_temperature_acceptance_ratio =
            max(max_temperature_acceptance_ratio, maximum(errors.temperature_acceptance_ratios))
        max_relative_internal_energy_error =
            max(max_relative_internal_energy_error, maximum(errors.relative_energy_errors))
        max_relative_cice_enthalpy_error =
            max(max_relative_cice_enthalpy_error, maximum(errors.relative_enthalpy_errors))
        max_relative_column_energy_error =
            max(max_relative_column_energy_error, errors.relative_column_energy_error)

        current_temperature = predicted_temperature
        current_record = target_record

        moving_temperature =
            fixed_grid_forced_step(moving_record, moving_temperature,
                                   target_record, relation, energy_transport;
                                   remap_to_target_thickness = true)

        moving_errors = replay_profile_errors(moving_temperature, target_record, relation, energy_floor)
        moving_final_abs_temperature_error = maximum(moving_errors.absolute_temperature_errors)
        moving_final_relative_temperature_error = maximum(moving_errors.relative_temperature_errors)
        moving_max_abs_temperature_error =
            max(moving_max_abs_temperature_error, moving_final_abs_temperature_error)
        moving_max_relative_temperature_error =
            max(moving_max_relative_temperature_error, moving_final_relative_temperature_error)
        moving_max_temperature_acceptance_ratio =
            max(moving_max_temperature_acceptance_ratio,
                maximum(moving_errors.temperature_acceptance_ratios))
        moving_max_relative_internal_energy_error =
            max(moving_max_relative_internal_energy_error, maximum(moving_errors.relative_energy_errors))
        moving_max_relative_cice_enthalpy_error =
            max(moving_max_relative_cice_enthalpy_error, maximum(moving_errors.relative_enthalpy_errors))
        moving_max_relative_column_energy_error =
            max(moving_max_relative_column_energy_error, moving_errors.relative_column_energy_error)
        moving_record = target_record

        moving_metric_temperature =
            moving_metric_forced_step(moving_metric_record,
                                      moving_metric_temperature,
                                      target_record,
                                      relation,
                                      energy_transport)

        moving_metric_errors = replay_profile_errors(moving_metric_temperature,
                                                     target_record,
                                                     relation,
                                                     energy_floor)
        moving_metric_final_abs_temperature_error =
            maximum(moving_metric_errors.absolute_temperature_errors)
        moving_metric_final_relative_temperature_error =
            maximum(moving_metric_errors.relative_temperature_errors)
        moving_metric_max_abs_temperature_error =
            max(moving_metric_max_abs_temperature_error,
                moving_metric_final_abs_temperature_error)
        moving_metric_max_relative_temperature_error =
            max(moving_metric_max_relative_temperature_error,
                moving_metric_final_relative_temperature_error)
        moving_metric_max_temperature_acceptance_ratio =
            max(moving_metric_max_temperature_acceptance_ratio,
                maximum(moving_metric_errors.temperature_acceptance_ratios))
        moving_metric_max_relative_internal_energy_error =
            max(moving_metric_max_relative_internal_energy_error,
                maximum(moving_metric_errors.relative_energy_errors))
        moving_metric_max_relative_cice_enthalpy_error =
            max(moving_metric_max_relative_cice_enthalpy_error,
                maximum(moving_metric_errors.relative_enthalpy_errors))
        moving_metric_max_relative_column_energy_error =
            max(moving_metric_max_relative_column_energy_error,
                moving_metric_errors.relative_column_energy_error)
        moving_metric_record = target_record

        prognostic_step =
            prognostic_thickness_forced_step(prognostic_record,
                                             prognostic_temperature,
                                             prognostic_thickness,
                                             target_record,
                                             relation,
                                             energy_transport;
                                             fixed_thickness = fixed_thickness_case)

        prognostic_temperature = prognostic_step.temperature
        prognostic_thickness = prognostic_step.thickness
        prognostic_top_melt += max(0.0, -prognostic_step.surface_thickness_change)
        prognostic_basal_growth += max(0.0, prognostic_step.basal_thickness_change)
        prognostic_bottom_melt += max(0.0, -prognostic_step.basal_thickness_change)

        prognostic_errors =
            replay_profile_errors(prognostic_temperature, target_record, relation, energy_floor)
        prognostic_max_temperature_acceptance_ratio =
            max(prognostic_max_temperature_acceptance_ratio,
                maximum(prognostic_errors.temperature_acceptance_ratios))
        prognostic_max_relative_cice_enthalpy_error =
            max(prognostic_max_relative_cice_enthalpy_error,
                maximum(prognostic_errors.relative_enthalpy_errors))
        prognostic_max_relative_column_energy_error =
            max(prognostic_max_relative_column_energy_error,
                prognostic_errors.relative_column_energy_error)
        prognostic_record = target_record

        icepack_matrix_prognostic_step =
            icepack_matrix_prognostic_thickness_forced_step(icepack_matrix_prognostic_record,
                                                            icepack_matrix_prognostic_temperature,
                                                              icepack_matrix_prognostic_thickness,
                                                              target_record,
                                                              relation,
                                                              energy_transport;
                                                              fixed_thickness = fixed_thickness_case)

        icepack_matrix_prognostic_temperature = icepack_matrix_prognostic_step.temperature
        icepack_matrix_prognostic_thickness = icepack_matrix_prognostic_step.thickness
        icepack_matrix_prognostic_top_melt +=
            max(0.0, -icepack_matrix_prognostic_step.surface_thickness_change)
        icepack_matrix_prognostic_basal_growth +=
            max(0.0, icepack_matrix_prognostic_step.basal_thickness_change)
        icepack_matrix_prognostic_bottom_melt +=
            max(0.0, -icepack_matrix_prognostic_step.basal_thickness_change)

        icepack_matrix_prognostic_errors =
            replay_profile_errors(icepack_matrix_prognostic_temperature,
                                  target_record,
                                  relation,
                                  energy_floor)
        icepack_matrix_prognostic_max_temperature_acceptance_ratio =
            max(icepack_matrix_prognostic_max_temperature_acceptance_ratio,
                maximum(icepack_matrix_prognostic_errors.temperature_acceptance_ratios))
        icepack_matrix_prognostic_max_relative_cice_enthalpy_error =
            max(icepack_matrix_prognostic_max_relative_cice_enthalpy_error,
                maximum(icepack_matrix_prognostic_errors.relative_enthalpy_errors))
        icepack_matrix_prognostic_max_relative_column_energy_error =
            max(icepack_matrix_prognostic_max_relative_column_energy_error,
                icepack_matrix_prognostic_errors.relative_column_energy_error)
        icepack_matrix_prognostic_record = target_record

        icepack_matrix_moving_temperature =
            icepack_temperature_matrix_replay_step(icepack_matrix_moving_record,
                                                   target_record,
                                                   relation,
                                                   energy_transport;
                                                   start_temperature = icepack_matrix_moving_temperature,
                                                   remap_to_target_thickness = true)

        icepack_matrix_moving_errors =
            replay_profile_errors(icepack_matrix_moving_temperature, target_record, relation, energy_floor)
        icepack_matrix_moving_final_abs_temperature_error =
            maximum(icepack_matrix_moving_errors.absolute_temperature_errors)
        icepack_matrix_moving_final_relative_temperature_error =
            maximum(icepack_matrix_moving_errors.relative_temperature_errors)
        icepack_matrix_moving_final_temperature_acceptance_ratio =
            maximum(icepack_matrix_moving_errors.temperature_acceptance_ratios)
        icepack_matrix_moving_final_top_temperature_error =
            first(icepack_matrix_moving_temperature) - first(target_record.Tice)
        icepack_matrix_moving_final_bottom_temperature_error =
            last(icepack_matrix_moving_temperature) - last(target_record.Tice)

        interval_max_abs_temperature_error, interval_max_error_layer =
            findmax(icepack_matrix_moving_errors.absolute_temperature_errors)
        if interval_max_abs_temperature_error > icepack_matrix_moving_max_abs_temperature_error
            icepack_matrix_moving_max_error_time_days = target_record.time_days
            icepack_matrix_moving_max_error_layer_from_top = interval_max_error_layer
            icepack_matrix_moving_max_error_predicted_temperature =
                icepack_matrix_moving_temperature[interval_max_error_layer]
            icepack_matrix_moving_max_error_cice_temperature =
                target_record.Tice[interval_max_error_layer]
        end

        icepack_matrix_moving_max_abs_temperature_error =
            max(icepack_matrix_moving_max_abs_temperature_error,
                icepack_matrix_moving_final_abs_temperature_error)
        icepack_matrix_moving_max_relative_temperature_error =
            max(icepack_matrix_moving_max_relative_temperature_error,
                icepack_matrix_moving_final_relative_temperature_error)
        icepack_matrix_moving_max_temperature_acceptance_ratio =
            max(icepack_matrix_moving_max_temperature_acceptance_ratio,
                icepack_matrix_moving_final_temperature_acceptance_ratio)
        icepack_matrix_moving_max_relative_internal_energy_error =
            max(icepack_matrix_moving_max_relative_internal_energy_error,
                maximum(icepack_matrix_moving_errors.relative_energy_errors))
        icepack_matrix_moving_max_relative_cice_enthalpy_error =
            max(icepack_matrix_moving_max_relative_cice_enthalpy_error,
                maximum(icepack_matrix_moving_errors.relative_enthalpy_errors))
        icepack_matrix_moving_max_relative_column_energy_error =
            max(icepack_matrix_moving_max_relative_column_energy_error,
                icepack_matrix_moving_errors.relative_column_energy_error)
        icepack_matrix_moving_record = target_record
    end

    interval_starts = records[1:end-1]
    interval_targets = records[2:end]

    for (n, (start_record, target_record)) in enumerate(zip(interval_starts, interval_targets))
        predicted_temperature =
            fixed_grid_forced_step(start_record, start_record.Tice,
                                   target_record, relation, energy_transport)

        reset_errors = replay_profile_errors(predicted_temperature, target_record, relation, energy_floor)
        interval_abs_error = maximum(reset_errors.absolute_temperature_errors)
        interval_relative_error = maximum(reset_errors.relative_temperature_errors)

        n == 1 && (reset_first_abs_temperature_error = interval_abs_error)
        n == 1 && (reset_first_relative_temperature_error = interval_relative_error)

        reset_max_abs_temperature_error =
            max(reset_max_abs_temperature_error, interval_abs_error)
        reset_max_relative_temperature_error =
            max(reset_max_relative_temperature_error, interval_relative_error)
        reset_max_temperature_acceptance_ratio =
            max(reset_max_temperature_acceptance_ratio,
                maximum(reset_errors.temperature_acceptance_ratios))
        reset_max_relative_internal_energy_error =
            max(reset_max_relative_internal_energy_error, maximum(reset_errors.relative_energy_errors))
        reset_max_relative_cice_enthalpy_error =
            max(reset_max_relative_cice_enthalpy_error, maximum(reset_errors.relative_enthalpy_errors))
        reset_max_relative_column_energy_error =
            max(reset_max_relative_column_energy_error, reset_errors.relative_column_energy_error)

        icepack_matrix_temperature =
            icepack_temperature_matrix_replay_step(start_record,
                                                   target_record,
                                                   relation,
                                                   energy_transport)
        icepack_matrix_errors =
            replay_profile_errors(icepack_matrix_temperature, target_record, relation, energy_floor)
        icepack_matrix_max_abs_temperature_error =
            max(icepack_matrix_max_abs_temperature_error,
                maximum(icepack_matrix_errors.absolute_temperature_errors))
        icepack_matrix_max_relative_temperature_error =
            max(icepack_matrix_max_relative_temperature_error,
                maximum(icepack_matrix_errors.relative_temperature_errors))
        icepack_matrix_max_temperature_acceptance_ratio =
            max(icepack_matrix_max_temperature_acceptance_ratio,
                maximum(icepack_matrix_errors.temperature_acceptance_ratios))
        icepack_matrix_max_relative_internal_energy_error =
            max(icepack_matrix_max_relative_internal_energy_error,
                maximum(icepack_matrix_errors.relative_energy_errors))
        icepack_matrix_max_relative_cice_enthalpy_error =
            max(icepack_matrix_max_relative_cice_enthalpy_error,
                maximum(icepack_matrix_errors.relative_enthalpy_errors))
        icepack_matrix_max_relative_column_energy_error =
            max(icepack_matrix_max_relative_column_energy_error,
                icepack_matrix_errors.relative_column_energy_error)

        if !isapprox(target_record.hi, start_record.hi; atol = 1e-12, rtol = 0)
            icepack_matrix_remapped_temperature =
                icepack_temperature_matrix_replay_step(start_record,
                                                       target_record,
                                                       relation,
                                                       energy_transport;
                                                       remap_to_target_thickness = true)
            icepack_matrix_remapped_errors =
                replay_profile_errors(icepack_matrix_remapped_temperature, target_record, relation, energy_floor)
            icepack_matrix_remapped_reset_max_abs_temperature_error =
                max(icepack_matrix_remapped_reset_max_abs_temperature_error,
                    maximum(icepack_matrix_remapped_errors.absolute_temperature_errors))
            icepack_matrix_remapped_reset_max_relative_temperature_error =
                max(icepack_matrix_remapped_reset_max_relative_temperature_error,
                    maximum(icepack_matrix_remapped_errors.relative_temperature_errors))
            icepack_matrix_remapped_reset_max_temperature_acceptance_ratio =
                max(icepack_matrix_remapped_reset_max_temperature_acceptance_ratio,
                    maximum(icepack_matrix_remapped_errors.temperature_acceptance_ratios))
            icepack_matrix_remapped_reset_max_relative_internal_energy_error =
                max(icepack_matrix_remapped_reset_max_relative_internal_energy_error,
                    maximum(icepack_matrix_remapped_errors.relative_energy_errors))
            icepack_matrix_remapped_reset_max_relative_cice_enthalpy_error =
                max(icepack_matrix_remapped_reset_max_relative_cice_enthalpy_error,
                    maximum(icepack_matrix_remapped_errors.relative_enthalpy_errors))
            icepack_matrix_remapped_reset_max_relative_column_energy_error =
                max(icepack_matrix_remapped_reset_max_relative_column_energy_error,
                    icepack_matrix_remapped_errors.relative_column_energy_error)
            icepack_matrix_remapped_reset_intervals += 1

            thickness_remapped_temperature =
                fixed_grid_forced_step(start_record, start_record.Tice,
                                       target_record, relation, energy_transport;
                                       remap_to_target_thickness = true)

            thickness_remapped_errors =
                replay_profile_errors(thickness_remapped_temperature, target_record, relation, energy_floor)
            thickness_remapped_reset_max_abs_temperature_error =
                max(thickness_remapped_reset_max_abs_temperature_error,
                    maximum(thickness_remapped_errors.absolute_temperature_errors))
            thickness_remapped_reset_max_relative_temperature_error =
                max(thickness_remapped_reset_max_relative_temperature_error,
                    maximum(thickness_remapped_errors.relative_temperature_errors))
            thickness_remapped_reset_max_temperature_acceptance_ratio =
                max(thickness_remapped_reset_max_temperature_acceptance_ratio,
                    maximum(thickness_remapped_errors.temperature_acceptance_ratios))
            thickness_remapped_reset_max_relative_internal_energy_error =
                max(thickness_remapped_reset_max_relative_internal_energy_error,
                    maximum(thickness_remapped_errors.relative_energy_errors))
            thickness_remapped_reset_max_relative_cice_enthalpy_error =
                max(thickness_remapped_reset_max_relative_cice_enthalpy_error,
                    maximum(thickness_remapped_errors.relative_enthalpy_errors))
            thickness_remapped_reset_max_relative_column_energy_error =
                max(thickness_remapped_reset_max_relative_column_energy_error,
                    thickness_remapped_errors.relative_column_energy_error)
            thickness_remapped_reset_intervals += 1
        end

        if target_record.hi < start_record.hi
            remapped_temperature =
                fixed_grid_forced_step(start_record, start_record.Tice,
                                       target_record, relation, energy_transport;
                                       remap_to_target_thickness = true)

            remapped_errors = replay_profile_errors(remapped_temperature, target_record, relation, energy_floor)
            remapped_reset_max_abs_temperature_error =
                max(remapped_reset_max_abs_temperature_error,
                    maximum(remapped_errors.absolute_temperature_errors))
            remapped_reset_max_relative_temperature_error =
                max(remapped_reset_max_relative_temperature_error,
                    maximum(remapped_errors.relative_temperature_errors))
            remapped_reset_max_temperature_acceptance_ratio =
                max(remapped_reset_max_temperature_acceptance_ratio,
                    maximum(remapped_errors.temperature_acceptance_ratios))
            remapped_reset_max_relative_internal_energy_error =
                max(remapped_reset_max_relative_internal_energy_error,
                    maximum(remapped_errors.relative_energy_errors))
            remapped_reset_max_relative_cice_enthalpy_error =
                max(remapped_reset_max_relative_cice_enthalpy_error,
                    maximum(remapped_errors.relative_enthalpy_errors))
            remapped_reset_max_relative_column_energy_error =
                max(remapped_reset_max_relative_column_energy_error,
                    remapped_errors.relative_column_energy_error)
            remapped_reset_intervals += 1
        end
    end

    remapped_reset_status =
        remapped_reset_intervals == 0 ? "not_run" :
        remapped_reset_max_relative_temperature_error <= 0.01 ? "pass" : "fail"

    thickness_remapped_reset_status =
        thickness_remapped_reset_intervals == 0 ? "not_run" :
        thickness_remapped_reset_max_relative_temperature_error <= 0.01 ? "pass" : "fail"

    icepack_matrix_remapped_reset_status =
        icepack_matrix_remapped_reset_intervals == 0 ? "not_run" :
        icepack_matrix_remapped_reset_max_relative_temperature_error <= 0.01 ? "pass" : "fail"

    thickness_rates = cice_thickness_rate_integrals(records)
    prognostic_final_hi_error = abs(prognostic_thickness - last(records).hi)
    prognostic_top_melt_error = abs(prognostic_top_melt - thickness_rates.integrated_top_melt)
    prognostic_basal_growth_error = abs(prognostic_basal_growth - thickness_rates.integrated_basal_growth)
    prognostic_bottom_melt_error = abs(prognostic_bottom_melt - thickness_rates.integrated_bottom_melt)
    prognostic_final_hi_acceptance_ratio =
        thickness_acceptance_ratio(prognostic_thickness, last(records).hi)
    prognostic_top_melt_acceptance_ratio =
        thickness_acceptance_ratio(prognostic_top_melt, thickness_rates.integrated_top_melt)
    prognostic_basal_growth_acceptance_ratio =
        thickness_acceptance_ratio(prognostic_basal_growth, thickness_rates.integrated_basal_growth)
    prognostic_bottom_melt_acceptance_ratio =
        thickness_acceptance_ratio(prognostic_bottom_melt, thickness_rates.integrated_bottom_melt)
    fixed_thickness_drift = CASE_ID == "cold_conductive_relaxation" ?
                            abs(prognostic_thickness - initial_record.hi) :
                            0.0
    fixed_thickness_acceptance_ratio =
        CASE_ID == "cold_conductive_relaxation" ? fixed_thickness_drift / 1e-6 : 0.0
    prognostic_thickness_acceptance_ratio =
        maximum((prognostic_final_hi_acceptance_ratio,
                 prognostic_top_melt_acceptance_ratio,
                 prognostic_basal_growth_acceptance_ratio,
                 prognostic_bottom_melt_acceptance_ratio,
                 fixed_thickness_acceptance_ratio))
    prognostic_thermal_acceptance_ratio =
        maximum((prognostic_max_temperature_acceptance_ratio,
                 prognostic_max_relative_cice_enthalpy_error / 0.01,
                 prognostic_max_relative_column_energy_error / 0.01))
    prognostic_acceptance_ratio =
        maximum((prognostic_thickness_acceptance_ratio,
                 prognostic_thermal_acceptance_ratio))

    icepack_matrix_prognostic_final_hi_error =
        abs(icepack_matrix_prognostic_thickness - last(records).hi)
    icepack_matrix_prognostic_top_melt_error =
        abs(icepack_matrix_prognostic_top_melt - thickness_rates.integrated_top_melt)
    icepack_matrix_prognostic_basal_growth_error =
        abs(icepack_matrix_prognostic_basal_growth - thickness_rates.integrated_basal_growth)
    icepack_matrix_prognostic_bottom_melt_error =
        abs(icepack_matrix_prognostic_bottom_melt - thickness_rates.integrated_bottom_melt)
    icepack_matrix_prognostic_final_hi_acceptance_ratio =
        thickness_acceptance_ratio(icepack_matrix_prognostic_thickness, last(records).hi)
    icepack_matrix_prognostic_top_melt_acceptance_ratio =
        thickness_acceptance_ratio(icepack_matrix_prognostic_top_melt, thickness_rates.integrated_top_melt)
    icepack_matrix_prognostic_basal_growth_acceptance_ratio =
        thickness_acceptance_ratio(icepack_matrix_prognostic_basal_growth, thickness_rates.integrated_basal_growth)
    icepack_matrix_prognostic_bottom_melt_acceptance_ratio =
        thickness_acceptance_ratio(icepack_matrix_prognostic_bottom_melt, thickness_rates.integrated_bottom_melt)
    icepack_matrix_prognostic_fixed_thickness_drift =
        CASE_ID == "cold_conductive_relaxation" ?
        abs(icepack_matrix_prognostic_thickness - initial_record.hi) :
        0.0
    icepack_matrix_prognostic_fixed_thickness_acceptance_ratio =
        CASE_ID == "cold_conductive_relaxation" ?
        icepack_matrix_prognostic_fixed_thickness_drift / 1e-6 :
        0.0
    icepack_matrix_prognostic_thickness_acceptance_ratio =
        maximum((icepack_matrix_prognostic_final_hi_acceptance_ratio,
                 icepack_matrix_prognostic_top_melt_acceptance_ratio,
                 icepack_matrix_prognostic_basal_growth_acceptance_ratio,
                 icepack_matrix_prognostic_bottom_melt_acceptance_ratio,
                 icepack_matrix_prognostic_fixed_thickness_acceptance_ratio))
    icepack_matrix_prognostic_thermal_acceptance_ratio =
        maximum((icepack_matrix_prognostic_max_temperature_acceptance_ratio,
                 icepack_matrix_prognostic_max_relative_cice_enthalpy_error / 0.01,
                 icepack_matrix_prognostic_max_relative_column_energy_error / 0.01))
    icepack_matrix_prognostic_acceptance_ratio =
        maximum((icepack_matrix_prognostic_thickness_acceptance_ratio,
                 icepack_matrix_prognostic_thermal_acceptance_ratio))

    return ReplayMetric[
        ReplayMetric(metric = "prognostic_thickness_forced_replay_acceptance_ratio",
                     status = prognostic_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Sequential forced replay that predicts thickness from the CICE-compatible top flux split and predicted basal conductive residual"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_thickness_acceptance_ratio",
                     status = prognostic_thickness_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_thickness_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Maximum acceptance ratio for final thickness, top melt, basal growth, bottom melt, and cold fixed-thickness drift"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_fixed_thickness_drift_m",
                     status = fixed_thickness_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = fixed_thickness_drift,
                     tolerance = CASE_ID == "cold_conductive_relaxation" ? 1e-6 : NaN,
                     evidence = "Fixed-thickness drift used for the cold conductive relaxation case"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_fixed_thickness_acceptance_ratio",
                     status = fixed_thickness_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = fixed_thickness_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Cold-case fixed-thickness acceptance ratio against the 1e-6 m requirement"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_temperature_acceptance_ratio",
                     status = prognostic_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate for the prognostic-thickness forced replay"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_max_relative_cice_enthalpy_error",
                     status = prognostic_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = prognostic_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy for the prognostic-thickness forced replay"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_max_relative_column_energy_error",
                     status = prognostic_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = prognostic_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Layer-integrated shifted BL99 internal energy for the prognostic-thickness forced replay"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_final_hi_m",
                     status = "diagnostic",
                     value = prognostic_thickness,
                     tolerance = NaN,
                     evidence = "Final ice thickness predicted by the prognostic-thickness forced replay"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_final_hi_error_m",
                     status = prognostic_final_hi_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_final_hi_error,
                     tolerance = last(records).hi < 0.1 ? 1e-3 : 0.01 * abs(last(records).hi),
                     evidence = "Absolute error relative to final CICE hi"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_final_hi_acceptance_ratio",
                     status = prognostic_final_hi_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_final_hi_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Final-thickness acceptance ratio against CICE hi"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_top_melt_m",
                     status = "diagnostic",
                     value = prognostic_top_melt,
                     tolerance = NaN,
                     evidence = "Cumulative top melt predicted from the top flux split fcondtop_ai - fsurf_ai"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_top_melt_error_m",
                     status = prognostic_top_melt_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_top_melt_error,
                     tolerance = thickness_rates.integrated_top_melt < 0.1 ? 1e-3 : 0.01 * abs(thickness_rates.integrated_top_melt),
                     evidence = "Absolute top-melt error relative to integrated CICE meltt"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_top_melt_acceptance_ratio",
                     status = prognostic_top_melt_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_top_melt_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Top-melt acceptance ratio against integrated CICE meltt"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_basal_growth_m",
                     status = "diagnostic",
                     value = prognostic_basal_growth,
                     tolerance = NaN,
                     evidence = "Cumulative basal growth predicted from the bottom conductive residual"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_basal_growth_error_m",
                     status = prognostic_basal_growth_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_basal_growth_error,
                     tolerance = thickness_rates.integrated_basal_growth < 0.1 ? 1e-3 : 0.01 * abs(thickness_rates.integrated_basal_growth),
                     evidence = "Absolute basal-growth error relative to integrated CICE congel"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_basal_growth_acceptance_ratio",
                     status = prognostic_basal_growth_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_basal_growth_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Basal-growth acceptance ratio against integrated CICE congel"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_bottom_melt_m",
                     status = "diagnostic",
                     value = prognostic_bottom_melt,
                     tolerance = NaN,
                     evidence = "Cumulative bottom melt predicted from the bottom conductive residual"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_bottom_melt_error_m",
                     status = prognostic_bottom_melt_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_bottom_melt_error,
                     tolerance = thickness_rates.integrated_bottom_melt < 0.1 ? 1e-3 : 0.01 * abs(thickness_rates.integrated_bottom_melt),
                     evidence = "Absolute bottom-melt error relative to integrated CICE meltb"),
        ReplayMetric(metric = "prognostic_thickness_forced_replay_bottom_melt_acceptance_ratio",
                     status = prognostic_bottom_melt_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = prognostic_bottom_melt_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Bottom-melt acceptance ratio against integrated CICE meltb"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_acceptance_ratio",
                     status = icepack_matrix_prognostic_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Sequential source-level BL99 matrix replay that predicts thickness from top flux split and basal conductive residual"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_thickness_acceptance_ratio",
                     status = icepack_matrix_prognostic_thickness_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_thickness_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Source-level maximum acceptance ratio for final thickness, melt/growth, and cold fixed-thickness drift"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_fixed_thickness_drift_m",
                     status = icepack_matrix_prognostic_fixed_thickness_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_fixed_thickness_drift,
                     tolerance = CASE_ID == "cold_conductive_relaxation" ? 1e-6 : NaN,
                     evidence = "Source-level fixed-thickness drift used for cold conductive relaxation"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_fixed_thickness_acceptance_ratio",
                     status = icepack_matrix_prognostic_fixed_thickness_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_fixed_thickness_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Source-level cold-case fixed-thickness acceptance ratio against the 1e-6 m requirement"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_temperature_acceptance_ratio",
                     status = icepack_matrix_prognostic_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate for the source-level prognostic-thickness forced replay"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_max_relative_cice_enthalpy_error",
                     status = icepack_matrix_prognostic_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy for the source-level prognostic-thickness forced replay"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_max_relative_column_energy_error",
                     status = icepack_matrix_prognostic_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Layer-integrated shifted BL99 internal energy for the source-level prognostic-thickness forced replay"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_final_hi_m",
                     status = "diagnostic",
                     value = icepack_matrix_prognostic_thickness,
                     tolerance = NaN,
                     evidence = "Final ice thickness predicted by the source-level prognostic-thickness forced replay"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_final_hi_error_m",
                     status = icepack_matrix_prognostic_final_hi_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_final_hi_error,
                     tolerance = last(records).hi < 0.1 ? 1e-3 : 0.01 * abs(last(records).hi),
                     evidence = "Source-level absolute error relative to final CICE hi"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_top_melt_m",
                     status = "diagnostic",
                     value = icepack_matrix_prognostic_top_melt,
                     tolerance = NaN,
                     evidence = "Source-level cumulative top melt predicted from the top flux split"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_top_melt_error_m",
                     status = icepack_matrix_prognostic_top_melt_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_top_melt_error,
                     tolerance = thickness_rates.integrated_top_melt < 0.1 ? 1e-3 : 0.01 * abs(thickness_rates.integrated_top_melt),
                     evidence = "Source-level absolute top-melt error relative to integrated CICE meltt"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_basal_growth_m",
                     status = "diagnostic",
                     value = icepack_matrix_prognostic_basal_growth,
                     tolerance = NaN,
                     evidence = "Source-level cumulative basal growth predicted from the bottom conductive residual"),
        ReplayMetric(metric = "icepack_matrix_prognostic_thickness_forced_replay_basal_growth_error_m",
                     status = icepack_matrix_prognostic_basal_growth_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_prognostic_basal_growth_error,
                     tolerance = thickness_rates.integrated_basal_growth < 0.1 ? 1e-3 : 0.01 * abs(thickness_rates.integrated_basal_growth),
                     evidence = "Source-level absolute basal-growth error relative to integrated CICE congel"),
        ReplayMetric(metric = "fixed_grid_forced_replay_max_abs_temperature_error_C",
                     status = "diagnostic",
                     value = max_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Hourly CICE fcondtop_ai and prescribed bottom ocean temperature, $(CONDUCTIVITY_VARIANT) conductivity, no snow or thickness prognostic"),
        ReplayMetric(metric = "fixed_grid_forced_replay_max_relative_temperature_error",
                     status = max_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = max_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Relative to abs(CICE T), floored by 1 C"),
        ReplayMetric(metric = "fixed_grid_forced_replay_temperature_acceptance_ratio",
                     status = max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate: 1% relative error or 0.02 C absolute fallback when |T_CICE| < 2 C"),
        ReplayMetric(metric = "fixed_grid_forced_replay_max_relative_internal_energy_error",
                     status = max_relative_internal_energy_error <= 0.01 ? "pass" : "fail",
                     value = max_relative_internal_energy_error,
                     tolerance = 0.01,
                     evidence = "Shifted BL99 internal energy from predicted T and target CICE S; relative floor is fixed from cold Case 1"),
        ReplayMetric(metric = "fixed_grid_forced_replay_max_relative_cice_enthalpy_error",
                     status = max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy q_i = E - rho_i L0, with floor fixed from cold Case 1"),
        ReplayMetric(metric = "fixed_grid_forced_replay_max_relative_column_energy_error",
                     status = max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Layer-integrated shifted BL99 internal energy over target CICE thickness"),
        ReplayMetric(metric = "fixed_grid_forced_replay_final_abs_temperature_error_C",
                     status = "diagnostic",
                     value = final_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Final hourly record at day 5"),
        ReplayMetric(metric = "fixed_grid_forced_replay_final_relative_temperature_error",
                     status = final_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = final_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Final hourly record at day 5"),
        ReplayMetric(metric = "cice_thickness_forced_replay_max_abs_temperature_error_C",
                     status = "diagnostic",
                     value = moving_max_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Hourly CICE flux forcing, CICE hi sequence, conservative energy remap, and BL99 basal-congelation fill enthalpy"),
        ReplayMetric(metric = "cice_thickness_forced_replay_max_relative_temperature_error",
                     status = moving_max_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = moving_max_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Relative to abs(CICE T), floored by 1 C; thickness is forced from CICE history"),
        ReplayMetric(metric = "cice_thickness_forced_replay_temperature_acceptance_ratio",
                     status = moving_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = moving_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate for CICE-thickness-forced remap"),
        ReplayMetric(metric = "cice_thickness_forced_replay_max_relative_internal_energy_error",
                     status = moving_max_relative_internal_energy_error <= 0.01 ? "pass" : "fail",
                     value = moving_max_relative_internal_energy_error,
                     tolerance = 0.01,
                     evidence = "Shifted BL99 internal energy after CICE-thickness-forced remap"),
        ReplayMetric(metric = "cice_thickness_forced_replay_max_relative_cice_enthalpy_error",
                     status = moving_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = moving_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy after CICE-thickness-forced remap"),
        ReplayMetric(metric = "cice_thickness_forced_replay_max_relative_column_energy_error",
                     status = moving_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = moving_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Layer-integrated shifted BL99 internal energy after CICE-thickness-forced remap"),
        ReplayMetric(metric = "cice_thickness_forced_replay_final_abs_temperature_error_C",
                     status = "diagnostic",
                     value = moving_final_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Final hourly record after CICE-thickness-forced remap"),
        ReplayMetric(metric = "cice_thickness_forced_replay_final_relative_temperature_error",
                     status = moving_final_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = moving_final_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Final hourly record after CICE-thickness-forced remap"),
        ReplayMetric(metric = "moving_metric_cice_thickness_forced_replay_max_abs_temperature_error_C",
                     status = "diagnostic",
                     value = moving_metric_max_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Hourly CICE flux forcing and CICE hi sequence integrated on MutableVerticalDiscretization with moving Jacobian metric"),
        ReplayMetric(metric = "moving_metric_cice_thickness_forced_replay_max_relative_temperature_error",
                     status = moving_metric_max_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = moving_metric_max_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Relative to abs(CICE T), floored by 1 C; CICE thickness is represented by the mutable vertical metric"),
        ReplayMetric(metric = "moving_metric_cice_thickness_forced_replay_temperature_acceptance_ratio",
                     status = moving_metric_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = moving_metric_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate for MutableVerticalDiscretization moving-metric replay"),
        ReplayMetric(metric = "moving_metric_cice_thickness_forced_replay_max_relative_internal_energy_error",
                     status = moving_metric_max_relative_internal_energy_error <= 0.01 ? "pass" : "fail",
                     value = moving_metric_max_relative_internal_energy_error,
                     tolerance = 0.01,
                     evidence = "Shifted BL99 internal energy from MutableVerticalDiscretization moving-metric replay"),
        ReplayMetric(metric = "moving_metric_cice_thickness_forced_replay_max_relative_cice_enthalpy_error",
                     status = moving_metric_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = moving_metric_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy from MutableVerticalDiscretization moving-metric replay"),
        ReplayMetric(metric = "moving_metric_cice_thickness_forced_replay_max_relative_column_energy_error",
                     status = moving_metric_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = moving_metric_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Layer-integrated shifted BL99 internal energy from MutableVerticalDiscretization moving-metric replay"),
        ReplayMetric(metric = "moving_metric_cice_thickness_forced_replay_final_abs_temperature_error_C",
                     status = "diagnostic",
                     value = moving_metric_final_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Final hourly record for MutableVerticalDiscretization moving-metric replay"),
        ReplayMetric(metric = "moving_metric_cice_thickness_forced_replay_final_relative_temperature_error",
                     status = moving_metric_final_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = moving_metric_final_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Final hourly record for MutableVerticalDiscretization moving-metric replay"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_max_abs_temperature_error_C",
                     status = "diagnostic",
                     value = icepack_matrix_moving_max_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Sequential source-level BL99 temperature matrix, CICE hi sequence, and conservative energy remap"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_temperature_error",
                     status = icepack_matrix_moving_max_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_moving_max_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Relative to abs(CICE T), floored by 1 C; source-level matrix with CICE thickness history"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_temperature_acceptance_ratio",
                     status = icepack_matrix_moving_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_moving_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate: 1% relative error or 0.02 C absolute fallback when |T_CICE| < 2 C"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_internal_energy_error",
                     status = icepack_matrix_moving_max_relative_internal_energy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_moving_max_relative_internal_energy_error,
                     tolerance = 0.01,
                     evidence = "Shifted BL99 internal energy for sequential source-level matrix with CICE thickness history"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_cice_enthalpy_error",
                     status = icepack_matrix_moving_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_moving_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy for sequential source-level matrix with CICE thickness history"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_column_energy_error",
                     status = icepack_matrix_moving_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_moving_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Layer-integrated shifted BL99 internal energy for sequential source-level matrix with CICE thickness history"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_final_abs_temperature_error_C",
                     status = "diagnostic",
                     value = icepack_matrix_moving_final_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Final hourly record for sequential source-level matrix with CICE thickness history"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_final_relative_temperature_error",
                     status = icepack_matrix_moving_final_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_moving_final_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Final hourly record for sequential source-level matrix with CICE thickness history"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_final_temperature_acceptance_ratio",
                     status = icepack_matrix_moving_final_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_moving_final_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Final hourly record under the official temperature gate with near-zero absolute fallback"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_max_error_time_days",
                     status = "diagnostic",
                     value = icepack_matrix_moving_max_error_time_days,
                     tolerance = NaN,
                     evidence = "CICE history time where the sequential source-level matrix reached its maximum absolute temperature error"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_max_error_layer_from_top",
                     status = "diagnostic",
                     value = icepack_matrix_moving_max_error_layer_from_top,
                     tolerance = NaN,
                     evidence = "One-based CICE ice layer index, top to bottom, for the maximum sequential source-level matrix temperature error"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_max_error_predicted_temperature_C",
                     status = "diagnostic",
                     value = icepack_matrix_moving_max_error_predicted_temperature,
                     tolerance = NaN,
                     evidence = "Predicted temperature at the maximum-error time and layer"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_max_error_cice_temperature_C",
                     status = "diagnostic",
                     value = icepack_matrix_moving_max_error_cice_temperature,
                     tolerance = NaN,
                     evidence = "CICE temperature at the maximum-error time and layer"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_final_top_temperature_error_C",
                     status = "diagnostic",
                     value = icepack_matrix_moving_final_top_temperature_error,
                     tolerance = NaN,
                     evidence = "Predicted minus CICE top-layer temperature at the final hourly record"),
        ReplayMetric(metric = "icepack_temperature_matrix_cice_thickness_forced_replay_final_bottom_temperature_error_C",
                     status = "diagnostic",
                     value = icepack_matrix_moving_final_bottom_temperature_error,
                     tolerance = NaN,
                     evidence = "Predicted minus CICE bottom-layer temperature at the final hourly record"),
        ReplayMetric(metric = "reset_forced_step_max_abs_temperature_error_C",
                     status = "diagnostic",
                     value = reset_max_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Each interval starts from parsed CICE state"),
        ReplayMetric(metric = "reset_forced_step_max_relative_temperature_error",
                     status = reset_max_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = reset_max_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Each interval starts from parsed CICE state"),
        ReplayMetric(metric = "reset_forced_step_temperature_acceptance_ratio",
                     status = reset_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = reset_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate for reset-each-interval fixed-grid replay"),
        ReplayMetric(metric = "reset_forced_step_max_relative_internal_energy_error",
                     status = reset_max_relative_internal_energy_error <= 0.01 ? "pass" : "fail",
                     value = reset_max_relative_internal_energy_error,
                     tolerance = 0.01,
                     evidence = "Each interval starts from parsed CICE state; shifted BL99 internal energy"),
        ReplayMetric(metric = "reset_forced_step_max_relative_cice_enthalpy_error",
                     status = reset_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = reset_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "Each interval starts from parsed CICE state; CICE negative enthalpy"),
        ReplayMetric(metric = "reset_forced_step_max_relative_column_energy_error",
                     status = reset_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = reset_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Each interval starts from parsed CICE state; target-thickness column integral"),
        ReplayMetric(metric = "reset_forced_step_first_abs_temperature_error_C",
                     status = "diagnostic",
                     value = reset_first_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Initial to first-hour interval"),
        ReplayMetric(metric = "reset_forced_step_first_relative_temperature_error",
                     status = reset_first_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = reset_first_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Initial to first-hour interval"),
        ReplayMetric(metric = "icepack_temperature_matrix_reset_max_abs_temperature_error_C",
                     status = "diagnostic",
                     value = icepack_matrix_max_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Source-level BL99 temperature matrix with CICE top conductive flux and bottom ocean temperature"),
        ReplayMetric(metric = "icepack_temperature_matrix_reset_max_relative_temperature_error",
                     status = icepack_matrix_max_relative_temperature_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_max_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Source-level BL99 temperature matrix reset replay"),
        ReplayMetric(metric = "icepack_temperature_matrix_reset_temperature_acceptance_ratio",
                     status = icepack_matrix_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate for source-level BL99 temperature matrix reset replay"),
        ReplayMetric(metric = "icepack_temperature_matrix_reset_max_relative_internal_energy_error",
                     status = icepack_matrix_max_relative_internal_energy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_max_relative_internal_energy_error,
                     tolerance = 0.01,
                     evidence = "Shifted BL99 internal energy from source-level temperature-matrix replay"),
        ReplayMetric(metric = "icepack_temperature_matrix_reset_max_relative_cice_enthalpy_error",
                     status = icepack_matrix_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy from source-level temperature-matrix replay"),
        ReplayMetric(metric = "icepack_temperature_matrix_reset_max_relative_column_energy_error",
                     status = icepack_matrix_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Layer-integrated shifted BL99 internal energy from source-level temperature-matrix replay"),
        ReplayMetric(metric = "icepack_temperature_matrix_thickness_remapped_reset_max_abs_temperature_error_C",
                     status = icepack_matrix_remapped_reset_intervals == 0 ? "not_run" : "diagnostic",
                     value = icepack_matrix_remapped_reset_intervals == 0 ? NaN : icepack_matrix_remapped_reset_max_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Source-level BL99 temperature matrix replay followed by conservative CICE-style layer enthalpy remapping"),
        ReplayMetric(metric = "icepack_temperature_matrix_thickness_remapped_reset_max_relative_temperature_error",
                     status = icepack_matrix_remapped_reset_status,
                     value = icepack_matrix_remapped_reset_intervals == 0 ? NaN : icepack_matrix_remapped_reset_max_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Relative to abs(CICE T), floored by 1 C, after source-level temperature step and thickness remapping"),
        ReplayMetric(metric = "icepack_temperature_matrix_thickness_remapped_reset_temperature_acceptance_ratio",
                     status = icepack_matrix_remapped_reset_intervals == 0 ? "not_run" :
                              icepack_matrix_remapped_reset_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = icepack_matrix_remapped_reset_intervals == 0 ? NaN : icepack_matrix_remapped_reset_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate after source-level temperature step and thickness remapping"),
        ReplayMetric(metric = "icepack_temperature_matrix_thickness_remapped_reset_max_relative_internal_energy_error",
                     status = icepack_matrix_remapped_reset_intervals == 0 ? "not_run" :
                              icepack_matrix_remapped_reset_max_relative_internal_energy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_remapped_reset_intervals == 0 ? NaN : icepack_matrix_remapped_reset_max_relative_internal_energy_error,
                     tolerance = 0.01,
                     evidence = "Shifted BL99 internal energy after source-level temperature step and conservative thickness remapping"),
        ReplayMetric(metric = "icepack_temperature_matrix_thickness_remapped_reset_max_relative_cice_enthalpy_error",
                     status = icepack_matrix_remapped_reset_intervals == 0 ? "not_run" :
                              icepack_matrix_remapped_reset_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_remapped_reset_intervals == 0 ? NaN : icepack_matrix_remapped_reset_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy after source-level temperature step and conservative thickness remapping"),
        ReplayMetric(metric = "icepack_temperature_matrix_thickness_remapped_reset_max_relative_column_energy_error",
                     status = icepack_matrix_remapped_reset_intervals == 0 ? "not_run" :
                              icepack_matrix_remapped_reset_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = icepack_matrix_remapped_reset_intervals == 0 ? NaN : icepack_matrix_remapped_reset_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Column-integrated shifted BL99 internal energy after source-level temperature step and thickness remapping"),
        ReplayMetric(metric = "thickness_remapped_reset_max_abs_temperature_error_C",
                     status = thickness_remapped_reset_intervals == 0 ? "not_run" : "diagnostic",
                     value = thickness_remapped_reset_intervals == 0 ? NaN : thickness_remapped_reset_max_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Thickness-change intervals remapped conservatively from start thickness to target CICE thickness"),
        ReplayMetric(metric = "thickness_remapped_reset_max_relative_temperature_error",
                     status = thickness_remapped_reset_status,
                     value = thickness_remapped_reset_intervals == 0 ? NaN : thickness_remapped_reset_max_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Relative to abs(CICE T), floored by 1 C, after conservative thickness remap"),
        ReplayMetric(metric = "thickness_remapped_reset_temperature_acceptance_ratio",
                     status = thickness_remapped_reset_intervals == 0 ? "not_run" :
                              thickness_remapped_reset_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = thickness_remapped_reset_intervals == 0 ? NaN : thickness_remapped_reset_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate after conservative thickness remap"),
        ReplayMetric(metric = "thickness_remapped_reset_max_relative_internal_energy_error",
                     status = thickness_remapped_reset_intervals == 0 ? "not_run" :
                              thickness_remapped_reset_max_relative_internal_energy_error <= 0.01 ? "pass" : "fail",
                     value = thickness_remapped_reset_intervals == 0 ? NaN : thickness_remapped_reset_max_relative_internal_energy_error,
                     tolerance = 0.01,
                     evidence = "Shifted BL99 internal energy after conservative thickness remap"),
        ReplayMetric(metric = "thickness_remapped_reset_max_relative_cice_enthalpy_error",
                     status = thickness_remapped_reset_intervals == 0 ? "not_run" :
                              thickness_remapped_reset_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = thickness_remapped_reset_intervals == 0 ? NaN : thickness_remapped_reset_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy after conservative thickness remap"),
        ReplayMetric(metric = "thickness_remapped_reset_max_relative_column_energy_error",
                     status = thickness_remapped_reset_intervals == 0 ? "not_run" :
                              thickness_remapped_reset_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = thickness_remapped_reset_intervals == 0 ? NaN : thickness_remapped_reset_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Layer-integrated shifted BL99 internal energy after conservative thickness remap"),
        ReplayMetric(metric = "top_ablation_remapped_reset_max_abs_temperature_error_C",
                     status = remapped_reset_intervals == 0 ? "not_run" : "diagnostic",
                     value = remapped_reset_intervals == 0 ? NaN : remapped_reset_max_abs_temperature_error,
                     tolerance = NaN,
                     evidence = "Top-ablation intervals remapped conservatively from start thickness to target thickness"),
        ReplayMetric(metric = "top_ablation_remapped_reset_max_relative_temperature_error",
                     status = remapped_reset_status,
                     value = remapped_reset_intervals == 0 ? NaN : remapped_reset_max_relative_temperature_error,
                     tolerance = 0.01,
                     evidence = "Relative to abs(CICE T), floored by 1 C, after conservative top-ablation remap"),
        ReplayMetric(metric = "top_ablation_remapped_reset_temperature_acceptance_ratio",
                     status = remapped_reset_intervals == 0 ? "not_run" :
                              remapped_reset_max_temperature_acceptance_ratio <= 1 ? "pass" : "fail",
                     value = remapped_reset_intervals == 0 ? NaN : remapped_reset_max_temperature_acceptance_ratio,
                     tolerance = 1,
                     evidence = "Official temperature gate after conservative top-ablation remap"),
        ReplayMetric(metric = "top_ablation_remapped_reset_max_relative_internal_energy_error",
                     status = remapped_reset_intervals == 0 ? "not_run" :
                              remapped_reset_max_relative_internal_energy_error <= 0.01 ? "pass" : "fail",
                     value = remapped_reset_intervals == 0 ? NaN : remapped_reset_max_relative_internal_energy_error,
                     tolerance = 0.01,
                     evidence = "Shifted BL99 internal energy after conservative top-ablation remap"),
        ReplayMetric(metric = "top_ablation_remapped_reset_max_relative_cice_enthalpy_error",
                     status = remapped_reset_intervals == 0 ? "not_run" :
                              remapped_reset_max_relative_cice_enthalpy_error <= 0.01 ? "pass" : "fail",
                     value = remapped_reset_intervals == 0 ? NaN : remapped_reset_max_relative_cice_enthalpy_error,
                     tolerance = 0.01,
                     evidence = "CICE negative enthalpy after conservative top-ablation remap"),
        ReplayMetric(metric = "top_ablation_remapped_reset_max_relative_column_energy_error",
                     status = remapped_reset_intervals == 0 ? "not_run" :
                              remapped_reset_max_relative_column_energy_error <= 0.01 ? "pass" : "fail",
                     value = remapped_reset_intervals == 0 ? NaN : remapped_reset_max_relative_column_energy_error,
                     tolerance = 0.01,
                     evidence = "Layer-integrated shifted BL99 internal energy after conservative top-ablation remap"),
        ReplayMetric(metric = "first_hour_bottom_conductive_flux_W_m2",
                     status = "diagnostic",
                     value = first_bottom_conductive_flux,
                     tolerance = NaN,
                     evidence = "BL99 bottom conductance from parsed state; CICE siflcondbot is sparse in hourly history"),
        ReplayMetric(metric = "cice_internal_energy_relative_floor_J_m3",
                     status = "diagnostic",
                     value = energy_floor,
                     tolerance = NaN,
                     evidence = "Frozen relative-error floor from the cold conductive MU71 CICE history"),
        ReplayMetric(metric = "cice_enthalpy_relative_floor_J_m3",
                     status = "diagnostic",
                     value = CICE_ENTHALPY_RELATIVE_FLOOR,
                     tolerance = NaN,
                     evidence = "Frozen CICE negative-enthalpy relative-error floor from the cold conductive MU71 CICE history"),
    ]
end

function write_replay_summary(path, metrics)
    mkpath(dirname(path))

    open(path, "w") do file
        println(file, "metric,status,value,tolerance,evidence")
        for metric in metrics
            println(file, join(csv_field.((metric.metric,
                                           metric.status,
                                           format_csv_value(metric.value),
                                           format_csv_value(metric.tolerance),
                                           metric.evidence)), ','))
        end
    end
end

function print_replay_summary(metrics, path)
    println()
    println("ClimaSeaIce state replay")
    println("========================")
    println("metrics_csv: $path")

    for metric in metrics
        @printf("%-62s %-5s value=% .8e tolerance=% .8e\n",
                metric.metric, metric.status, metric.value, metric.tolerance)
    end
end

function print_forced_replay_summary(metrics, path)
    println()
    println("Fixed-grid forced replay")
    println("========================")
    println("metrics_csv: $path")

    for metric in metrics
        @printf("%-62s %-10s value=% .8e tolerance=% .8e\n",
                metric.metric, metric.status, metric.value, metric.tolerance)
    end
end

function metric_value(metrics, name)
    metric = metrics[findfirst(metric -> metric.metric == name, metrics)]
    return metric.value
end

function metric_status(metrics, name)
    metric = metrics[findfirst(metric -> metric.metric == name, metrics)]
    return metric.status
end

function read_csv_rows(path, header)
    isfile(path) || return Dict{String, String}[]
    lines = readlines(path)
    isempty(lines) && return Dict{String, String}[]
    file_header = split(first(lines), ',')
    rows = Dict{String, String}[]

    for line in lines[2:end]
        isempty(line) && continue
        file_row = Dict(file_header .=> split(line, ','; limit=length(file_header)))
        push!(rows, Dict(column => get(file_row, column, "") for column in header))
    end

    return rows
end

function write_csv_rows(path, header, rows)
    mkpath(dirname(path))
    open(path, "w") do file
        println(file, join(header, ','))
        for row in rows
            println(file, join([row[column] for column in header], ','))
        end
    end
end

function strict_case_comparison_status(path)
    header = case_comparison_header()
    rows = read_csv_rows(path, header)
    required_cases = Set(["cold_conductive_relaxation",
                          "surface_warming",
                          "surface_ablation",
                          "basal_growth"])
    required_conductivities = Set(["MU71", "bubbly"])
    observed = Set((row["case_id"], row["conduct"]) for row in rows)
    required = Set((case_id, conduct)
                   for case_id in required_cases
                   for conduct in required_conductivities)

    if !issubset(required, observed)
        return "incomplete"
    end

    return all(row["strict_bl99_validation_status"] == "pass"
               for row in rows
               if (row["case_id"], row["conduct"]) in required) ? "pass" : "incomplete"
end

case_comparison_header() =
    ["case_id",
     "conduct",
     "history_dir",
     "records",
     "hourly_records",
     "final_day",
     "initial_hi_m",
     "final_hi_m",
     "delta_hi_m",
     "initial_aice",
     "final_aice",
     "delta_aice",
     "initial_thermodynamic_hi_m",
     "final_thermodynamic_hi_m",
     "delta_thermodynamic_hi_m",
	     "thickness_loss_m",
	     "thickness_gain_m",
	     "initial_mean_Tice_C",
     "final_mean_Tice_C",
     "delta_mean_Tice_C",
     "max_congel_cm_per_day",
     "max_meltt_cm_per_day",
     "max_meltb_cm_per_day",
     "integrated_congel_m",
     "integrated_meltt_m",
     "integrated_meltb_m",
     "rate_integrated_delta_hi_m",
     "thickness_rate_budget_residual_m",
     "thickness_rate_budget_status",
     "case_transition_status",
     "prognostic_thickness_replay_status",
     "strict_bl99_validation_status",
     "prognostic_thickness_forced_replay_acceptance_ratio",
     "prognostic_thickness_forced_replay_thickness_acceptance_ratio",
     "prognostic_thickness_forced_replay_fixed_thickness_drift_m",
     "prognostic_thickness_forced_replay_fixed_thickness_acceptance_ratio",
     "prognostic_thickness_forced_replay_temperature_acceptance_ratio",
     "prognostic_thickness_forced_replay_max_relative_cice_enthalpy_error",
     "prognostic_thickness_forced_replay_max_relative_column_energy_error",
     "prognostic_thickness_forced_replay_final_hi_m",
     "prognostic_thickness_forced_replay_final_hi_error_m",
     "prognostic_thickness_forced_replay_final_hi_acceptance_ratio",
     "prognostic_thickness_forced_replay_top_melt_m",
     "prognostic_thickness_forced_replay_top_melt_error_m",
     "prognostic_thickness_forced_replay_top_melt_acceptance_ratio",
     "prognostic_thickness_forced_replay_basal_growth_m",
     "prognostic_thickness_forced_replay_basal_growth_error_m",
     "prognostic_thickness_forced_replay_basal_growth_acceptance_ratio",
     "prognostic_thickness_forced_replay_bottom_melt_m",
     "prognostic_thickness_forced_replay_bottom_melt_error_m",
     "prognostic_thickness_forced_replay_bottom_melt_acceptance_ratio",
     "icepack_matrix_prognostic_thickness_forced_status",
     "icepack_matrix_prognostic_thickness_forced_replay_acceptance_ratio",
     "icepack_matrix_prognostic_thickness_forced_replay_thickness_acceptance_ratio",
     "icepack_matrix_prognostic_thickness_forced_replay_fixed_thickness_drift_m",
     "icepack_matrix_prognostic_thickness_forced_replay_fixed_thickness_acceptance_ratio",
     "icepack_matrix_prognostic_thickness_forced_replay_temperature_acceptance_ratio",
     "icepack_matrix_prognostic_thickness_forced_replay_max_relative_cice_enthalpy_error",
     "icepack_matrix_prognostic_thickness_forced_replay_max_relative_column_energy_error",
     "icepack_matrix_prognostic_thickness_forced_replay_final_hi_m",
     "icepack_matrix_prognostic_thickness_forced_replay_final_hi_error_m",
     "icepack_matrix_prognostic_thickness_forced_replay_top_melt_m",
     "icepack_matrix_prognostic_thickness_forced_replay_top_melt_error_m",
     "icepack_matrix_prognostic_thickness_forced_replay_basal_growth_m",
     "icepack_matrix_prognostic_thickness_forced_replay_basal_growth_error_m",
     "state_replay_status",
     "reset_forced_step_status",
     "reset_forced_step_max_relative_temperature_error",
     "reset_forced_step_temperature_acceptance_status",
     "reset_forced_step_temperature_acceptance_ratio",
     "reset_forced_step_max_relative_internal_energy_error",
     "reset_forced_step_max_relative_cice_enthalpy_error",
     "reset_forced_step_max_relative_column_energy_error",
     "icepack_temperature_matrix_reset_status",
     "icepack_temperature_matrix_reset_max_relative_temperature_error",
     "icepack_temperature_matrix_reset_temperature_acceptance_status",
     "icepack_temperature_matrix_reset_temperature_acceptance_ratio",
     "icepack_temperature_matrix_reset_max_relative_internal_energy_error",
     "icepack_temperature_matrix_reset_max_relative_cice_enthalpy_error",
     "icepack_temperature_matrix_reset_max_relative_column_energy_error",
     "icepack_temperature_matrix_thickness_remapped_reset_status",
     "icepack_temperature_matrix_thickness_remapped_reset_max_relative_temperature_error",
     "icepack_temperature_matrix_thickness_remapped_reset_temperature_acceptance_status",
     "icepack_temperature_matrix_thickness_remapped_reset_temperature_acceptance_ratio",
     "icepack_temperature_matrix_thickness_remapped_reset_max_relative_internal_energy_error",
     "icepack_temperature_matrix_thickness_remapped_reset_max_relative_cice_enthalpy_error",
     "icepack_temperature_matrix_thickness_remapped_reset_max_relative_column_energy_error",
     "thickness_remapped_reset_status",
     "thickness_remapped_reset_max_relative_temperature_error",
     "thickness_remapped_reset_temperature_acceptance_status",
     "thickness_remapped_reset_temperature_acceptance_ratio",
     "thickness_remapped_reset_max_relative_internal_energy_error",
     "thickness_remapped_reset_max_relative_cice_enthalpy_error",
     "thickness_remapped_reset_max_relative_column_energy_error",
     "cice_thickness_forced_status",
     "cice_thickness_forced_max_relative_temperature_error",
     "cice_thickness_forced_temperature_acceptance_status",
     "cice_thickness_forced_temperature_acceptance_ratio",
     "cice_thickness_forced_max_relative_internal_energy_error",
     "cice_thickness_forced_max_relative_cice_enthalpy_error",
     "cice_thickness_forced_max_relative_column_energy_error",
     "moving_metric_cice_thickness_forced_status",
     "moving_metric_cice_thickness_forced_max_relative_temperature_error",
     "moving_metric_cice_thickness_forced_temperature_acceptance_status",
     "moving_metric_cice_thickness_forced_temperature_acceptance_ratio",
     "moving_metric_cice_thickness_forced_max_relative_internal_energy_error",
     "moving_metric_cice_thickness_forced_max_relative_cice_enthalpy_error",
     "moving_metric_cice_thickness_forced_max_relative_column_energy_error",
     "icepack_temperature_matrix_cice_thickness_forced_status",
     "icepack_temperature_matrix_cice_thickness_forced_max_relative_temperature_error",
     "icepack_temperature_matrix_cice_thickness_forced_temperature_acceptance_status",
     "icepack_temperature_matrix_cice_thickness_forced_temperature_acceptance_ratio",
     "icepack_temperature_matrix_cice_thickness_forced_max_relative_internal_energy_error",
     "icepack_temperature_matrix_cice_thickness_forced_max_relative_cice_enthalpy_error",
     "icepack_temperature_matrix_cice_thickness_forced_max_relative_column_energy_error",
     "top_ablation_remapped_reset_status",
     "top_ablation_remapped_reset_max_relative_temperature_error",
     "top_ablation_remapped_reset_temperature_acceptance_status",
     "top_ablation_remapped_reset_temperature_acceptance_ratio",
     "top_ablation_remapped_reset_max_relative_internal_energy_error",
     "top_ablation_remapped_reset_max_relative_cice_enthalpy_error",
     "top_ablation_remapped_reset_max_relative_column_energy_error",
     "free_running_status",
     "free_running_max_relative_temperature_error",
     "free_running_temperature_acceptance_status",
     "free_running_temperature_acceptance_ratio",
     "free_running_max_relative_internal_energy_error",
     "free_running_max_relative_cice_enthalpy_error",
     "free_running_max_relative_column_energy_error",
     "cice_summary_csv",
     "state_replay_csv",
     "forced_replay_csv"]

function case_transition_status(case_id, records)
    first_record = first(records)
    final_record = last(records)
    delta_hi = final_record.hi - first_record.hi
    thickness_loss = max(0.0, -delta_hi)
    max_growth = finite_maximum([record.congel for record in records])
    max_top_melt = finite_maximum([record.meltt for record in records])
    max_bottom_melt = finite_maximum([record.meltb for record in records])
    delta_mean_temperature = finite_mean(final_record.Tice) - finite_mean(first_record.Tice)

    if case_id == "cold_conductive_relaxation"
        no_growth_or_melt = max_growth <= 1e-8 && max_top_melt <= 1e-8 && max_bottom_melt <= 1e-8
        fixed_thickness = abs(delta_hi) <= 1e-5
        return no_growth_or_melt && fixed_thickness ? "pass" : "fail"
    elseif case_id == "surface_warming"
        return delta_mean_temperature > 0.0 && thickness_loss < 0.005 ? "pass" : "fail"
    elseif case_id == "surface_ablation"
        return thickness_loss >= 0.05 && max_top_melt > 0.0 ? "pass" : "fail"
    elseif case_id == "basal_growth"
        return delta_hi >= 0.01 && max_growth > 0.0 ? "pass" : "fail"
    else
        return "diagnostic"
    end
end

function upsert_case_comparison_summary!(path, history_dir, records, replay_metrics, forced_replay_metrics)
    header = case_comparison_header()
    case_id = cice_case_id(history_dir)
    conduct = CONDUCTIVITY_VARIANT
    rows = read_csv_rows(path, header)
    first_record = first(records)
    final_record = last(records)
    delta_hi = final_record.hi - first_record.hi
    initial_thermodynamic_hi = thermodynamic_ice_thickness(first_record)
    final_thermodynamic_hi = thermodynamic_ice_thickness(final_record)
    initial_mean_temperature = finite_mean(first_record.Tice)
    final_mean_temperature = finite_mean(final_record.Tice)
    thickness_rates = cice_thickness_rate_integrals(records)
    row = Dict(
        "case_id" => case_id,
        "conduct" => conduct,
        "history_dir" => history_dir,
        "records" => string(length(records)),
        "hourly_records" => string(count(record -> record.frequency == "hourly", records)),
        "final_day" => format_csv_value(final_record.time_days),
        "initial_hi_m" => format_csv_value(first_record.hi),
        "final_hi_m" => format_csv_value(final_record.hi),
        "delta_hi_m" => format_csv_value(delta_hi),
        "initial_aice" => format_csv_value(first_record.aice),
        "final_aice" => format_csv_value(final_record.aice),
        "delta_aice" => format_csv_value(final_record.aice - first_record.aice),
        "initial_thermodynamic_hi_m" => format_csv_value(initial_thermodynamic_hi),
        "final_thermodynamic_hi_m" => format_csv_value(final_thermodynamic_hi),
        "delta_thermodynamic_hi_m" => format_csv_value(final_thermodynamic_hi - initial_thermodynamic_hi),
        "thickness_loss_m" => format_csv_value(max(0.0, -delta_hi)),
        "thickness_gain_m" => format_csv_value(max(0.0, delta_hi)),
        "initial_mean_Tice_C" => format_csv_value(initial_mean_temperature),
        "final_mean_Tice_C" => format_csv_value(final_mean_temperature),
        "delta_mean_Tice_C" => format_csv_value(final_mean_temperature - initial_mean_temperature),
        "max_congel_cm_per_day" => format_csv_value(finite_maximum([record.congel for record in records])),
        "max_meltt_cm_per_day" => format_csv_value(finite_maximum([record.meltt for record in records])),
        "max_meltb_cm_per_day" => format_csv_value(finite_maximum([record.meltb for record in records])),
        "integrated_congel_m" => format_csv_value(thickness_rates.integrated_basal_growth),
        "integrated_meltt_m" => format_csv_value(thickness_rates.integrated_top_melt),
        "integrated_meltb_m" => format_csv_value(thickness_rates.integrated_bottom_melt),
        "rate_integrated_delta_hi_m" => format_csv_value(thickness_rates.rate_integrated_delta),
        "thickness_rate_budget_residual_m" => format_csv_value(thickness_rates.residual),
        "thickness_rate_budget_status" => cice_thickness_budget_status(records),
        "case_transition_status" => case_transition_status(case_id, records),
        "prognostic_thickness_replay_status" => metric_status(forced_replay_metrics, "prognostic_thickness_forced_replay_acceptance_ratio"),
        "strict_bl99_validation_status" => strict_bl99_validation_status(case_id, forced_replay_metrics),
        "prognostic_thickness_forced_replay_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_acceptance_ratio")),
        "prognostic_thickness_forced_replay_thickness_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_thickness_acceptance_ratio")),
        "prognostic_thickness_forced_replay_fixed_thickness_drift_m" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_fixed_thickness_drift_m")),
        "prognostic_thickness_forced_replay_fixed_thickness_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_fixed_thickness_acceptance_ratio")),
        "prognostic_thickness_forced_replay_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_temperature_acceptance_ratio")),
        "prognostic_thickness_forced_replay_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_max_relative_cice_enthalpy_error")),
        "prognostic_thickness_forced_replay_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_max_relative_column_energy_error")),
        "prognostic_thickness_forced_replay_final_hi_m" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_final_hi_m")),
        "prognostic_thickness_forced_replay_final_hi_error_m" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_final_hi_error_m")),
        "prognostic_thickness_forced_replay_final_hi_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_final_hi_acceptance_ratio")),
        "prognostic_thickness_forced_replay_top_melt_m" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_top_melt_m")),
        "prognostic_thickness_forced_replay_top_melt_error_m" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_top_melt_error_m")),
        "prognostic_thickness_forced_replay_top_melt_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_top_melt_acceptance_ratio")),
        "prognostic_thickness_forced_replay_basal_growth_m" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_basal_growth_m")),
        "prognostic_thickness_forced_replay_basal_growth_error_m" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_basal_growth_error_m")),
        "prognostic_thickness_forced_replay_basal_growth_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_basal_growth_acceptance_ratio")),
        "prognostic_thickness_forced_replay_bottom_melt_m" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_bottom_melt_m")),
        "prognostic_thickness_forced_replay_bottom_melt_error_m" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_bottom_melt_error_m")),
        "prognostic_thickness_forced_replay_bottom_melt_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "prognostic_thickness_forced_replay_bottom_melt_acceptance_ratio")),
        "icepack_matrix_prognostic_thickness_forced_status" => metric_status(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_acceptance_ratio"),
        "icepack_matrix_prognostic_thickness_forced_replay_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_acceptance_ratio")),
        "icepack_matrix_prognostic_thickness_forced_replay_thickness_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_thickness_acceptance_ratio")),
        "icepack_matrix_prognostic_thickness_forced_replay_fixed_thickness_drift_m" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_fixed_thickness_drift_m")),
        "icepack_matrix_prognostic_thickness_forced_replay_fixed_thickness_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_fixed_thickness_acceptance_ratio")),
        "icepack_matrix_prognostic_thickness_forced_replay_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_temperature_acceptance_ratio")),
        "icepack_matrix_prognostic_thickness_forced_replay_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_max_relative_cice_enthalpy_error")),
        "icepack_matrix_prognostic_thickness_forced_replay_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_max_relative_column_energy_error")),
        "icepack_matrix_prognostic_thickness_forced_replay_final_hi_m" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_final_hi_m")),
        "icepack_matrix_prognostic_thickness_forced_replay_final_hi_error_m" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_final_hi_error_m")),
        "icepack_matrix_prognostic_thickness_forced_replay_top_melt_m" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_top_melt_m")),
        "icepack_matrix_prognostic_thickness_forced_replay_top_melt_error_m" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_top_melt_error_m")),
        "icepack_matrix_prognostic_thickness_forced_replay_basal_growth_m" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_basal_growth_m")),
        "icepack_matrix_prognostic_thickness_forced_replay_basal_growth_error_m" => format_csv_value(metric_value(forced_replay_metrics, "icepack_matrix_prognostic_thickness_forced_replay_basal_growth_error_m")),
        "state_replay_status" => all(metric.status == "pass" for metric in replay_metrics) ? "pass" : "fail",
        "reset_forced_step_status" => metric_status(forced_replay_metrics, "reset_forced_step_max_relative_temperature_error"),
        "reset_forced_step_max_relative_temperature_error" => format_csv_value(metric_value(forced_replay_metrics, "reset_forced_step_max_relative_temperature_error")),
        "reset_forced_step_temperature_acceptance_status" => metric_status(forced_replay_metrics, "reset_forced_step_temperature_acceptance_ratio"),
        "reset_forced_step_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "reset_forced_step_temperature_acceptance_ratio")),
        "reset_forced_step_max_relative_internal_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "reset_forced_step_max_relative_internal_energy_error")),
        "reset_forced_step_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "reset_forced_step_max_relative_cice_enthalpy_error")),
        "reset_forced_step_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "reset_forced_step_max_relative_column_energy_error")),
        "icepack_temperature_matrix_reset_status" => metric_status(forced_replay_metrics, "icepack_temperature_matrix_reset_max_relative_temperature_error"),
        "icepack_temperature_matrix_reset_max_relative_temperature_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_reset_max_relative_temperature_error")),
        "icepack_temperature_matrix_reset_temperature_acceptance_status" => metric_status(forced_replay_metrics, "icepack_temperature_matrix_reset_temperature_acceptance_ratio"),
        "icepack_temperature_matrix_reset_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_reset_temperature_acceptance_ratio")),
        "icepack_temperature_matrix_reset_max_relative_internal_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_reset_max_relative_internal_energy_error")),
        "icepack_temperature_matrix_reset_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_reset_max_relative_cice_enthalpy_error")),
        "icepack_temperature_matrix_reset_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_reset_max_relative_column_energy_error")),
        "icepack_temperature_matrix_thickness_remapped_reset_status" => metric_status(forced_replay_metrics, "icepack_temperature_matrix_thickness_remapped_reset_max_relative_temperature_error"),
        "icepack_temperature_matrix_thickness_remapped_reset_max_relative_temperature_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_thickness_remapped_reset_max_relative_temperature_error")),
        "icepack_temperature_matrix_thickness_remapped_reset_temperature_acceptance_status" => metric_status(forced_replay_metrics, "icepack_temperature_matrix_thickness_remapped_reset_temperature_acceptance_ratio"),
        "icepack_temperature_matrix_thickness_remapped_reset_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_thickness_remapped_reset_temperature_acceptance_ratio")),
        "icepack_temperature_matrix_thickness_remapped_reset_max_relative_internal_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_thickness_remapped_reset_max_relative_internal_energy_error")),
        "icepack_temperature_matrix_thickness_remapped_reset_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_thickness_remapped_reset_max_relative_cice_enthalpy_error")),
        "icepack_temperature_matrix_thickness_remapped_reset_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_thickness_remapped_reset_max_relative_column_energy_error")),
        "thickness_remapped_reset_status" => metric_status(forced_replay_metrics, "thickness_remapped_reset_max_relative_temperature_error"),
        "thickness_remapped_reset_max_relative_temperature_error" => format_csv_value(metric_value(forced_replay_metrics, "thickness_remapped_reset_max_relative_temperature_error")),
        "thickness_remapped_reset_temperature_acceptance_status" => metric_status(forced_replay_metrics, "thickness_remapped_reset_temperature_acceptance_ratio"),
        "thickness_remapped_reset_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "thickness_remapped_reset_temperature_acceptance_ratio")),
        "thickness_remapped_reset_max_relative_internal_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "thickness_remapped_reset_max_relative_internal_energy_error")),
        "thickness_remapped_reset_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "thickness_remapped_reset_max_relative_cice_enthalpy_error")),
        "thickness_remapped_reset_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "thickness_remapped_reset_max_relative_column_energy_error")),
        "cice_thickness_forced_status" => metric_status(forced_replay_metrics, "cice_thickness_forced_replay_max_relative_temperature_error"),
        "cice_thickness_forced_max_relative_temperature_error" => format_csv_value(metric_value(forced_replay_metrics, "cice_thickness_forced_replay_max_relative_temperature_error")),
        "cice_thickness_forced_temperature_acceptance_status" => metric_status(forced_replay_metrics, "cice_thickness_forced_replay_temperature_acceptance_ratio"),
        "cice_thickness_forced_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "cice_thickness_forced_replay_temperature_acceptance_ratio")),
        "cice_thickness_forced_max_relative_internal_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "cice_thickness_forced_replay_max_relative_internal_energy_error")),
        "cice_thickness_forced_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "cice_thickness_forced_replay_max_relative_cice_enthalpy_error")),
        "cice_thickness_forced_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "cice_thickness_forced_replay_max_relative_column_energy_error")),
        "moving_metric_cice_thickness_forced_status" => metric_status(forced_replay_metrics, "moving_metric_cice_thickness_forced_replay_max_relative_temperature_error"),
        "moving_metric_cice_thickness_forced_max_relative_temperature_error" => format_csv_value(metric_value(forced_replay_metrics, "moving_metric_cice_thickness_forced_replay_max_relative_temperature_error")),
        "moving_metric_cice_thickness_forced_temperature_acceptance_status" => metric_status(forced_replay_metrics, "moving_metric_cice_thickness_forced_replay_temperature_acceptance_ratio"),
        "moving_metric_cice_thickness_forced_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "moving_metric_cice_thickness_forced_replay_temperature_acceptance_ratio")),
        "moving_metric_cice_thickness_forced_max_relative_internal_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "moving_metric_cice_thickness_forced_replay_max_relative_internal_energy_error")),
        "moving_metric_cice_thickness_forced_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "moving_metric_cice_thickness_forced_replay_max_relative_cice_enthalpy_error")),
        "moving_metric_cice_thickness_forced_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "moving_metric_cice_thickness_forced_replay_max_relative_column_energy_error")),
        "icepack_temperature_matrix_cice_thickness_forced_status" => metric_status(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_temperature_error"),
        "icepack_temperature_matrix_cice_thickness_forced_max_relative_temperature_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_temperature_error")),
        "icepack_temperature_matrix_cice_thickness_forced_temperature_acceptance_status" => metric_status(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_temperature_acceptance_ratio"),
        "icepack_temperature_matrix_cice_thickness_forced_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_temperature_acceptance_ratio")),
        "icepack_temperature_matrix_cice_thickness_forced_max_relative_internal_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_internal_energy_error")),
        "icepack_temperature_matrix_cice_thickness_forced_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_cice_enthalpy_error")),
        "icepack_temperature_matrix_cice_thickness_forced_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "icepack_temperature_matrix_cice_thickness_forced_replay_max_relative_column_energy_error")),
        "top_ablation_remapped_reset_status" => metric_status(forced_replay_metrics, "top_ablation_remapped_reset_max_relative_temperature_error"),
        "top_ablation_remapped_reset_max_relative_temperature_error" => format_csv_value(metric_value(forced_replay_metrics, "top_ablation_remapped_reset_max_relative_temperature_error")),
        "top_ablation_remapped_reset_temperature_acceptance_status" => metric_status(forced_replay_metrics, "top_ablation_remapped_reset_temperature_acceptance_ratio"),
        "top_ablation_remapped_reset_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "top_ablation_remapped_reset_temperature_acceptance_ratio")),
        "top_ablation_remapped_reset_max_relative_internal_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "top_ablation_remapped_reset_max_relative_internal_energy_error")),
        "top_ablation_remapped_reset_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "top_ablation_remapped_reset_max_relative_cice_enthalpy_error")),
        "top_ablation_remapped_reset_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "top_ablation_remapped_reset_max_relative_column_energy_error")),
        "free_running_status" => metric_status(forced_replay_metrics, "fixed_grid_forced_replay_max_relative_temperature_error"),
        "free_running_max_relative_temperature_error" => format_csv_value(metric_value(forced_replay_metrics, "fixed_grid_forced_replay_max_relative_temperature_error")),
        "free_running_temperature_acceptance_status" => metric_status(forced_replay_metrics, "fixed_grid_forced_replay_temperature_acceptance_ratio"),
        "free_running_temperature_acceptance_ratio" => format_csv_value(metric_value(forced_replay_metrics, "fixed_grid_forced_replay_temperature_acceptance_ratio")),
        "free_running_max_relative_internal_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "fixed_grid_forced_replay_max_relative_internal_energy_error")),
        "free_running_max_relative_cice_enthalpy_error" => format_csv_value(metric_value(forced_replay_metrics, "fixed_grid_forced_replay_max_relative_cice_enthalpy_error")),
        "free_running_max_relative_column_energy_error" => format_csv_value(metric_value(forced_replay_metrics, "fixed_grid_forced_replay_max_relative_column_energy_error")),
        "cice_summary_csv" => CICE_SUMMARY,
        "state_replay_csv" => CLIMASEAICE_REPLAY_SUMMARY,
        "forced_replay_csv" => FORCED_REPLAY_SUMMARY)

    index = findfirst(existing -> existing["case_id"] == case_id &&
                                  existing["conduct"] == conduct,
                      rows)
    if isnothing(index)
        push!(rows, row)
    else
        rows[index] = row
    end

    sort!(rows; by = row -> (row["case_id"], row["conduct"]))
    write_csv_rows(path, header, rows)
    return path
end

function print_current_case_status(history_dir)
    case_id = cice_case_id(history_dir)
    rows = read_csv_rows(CASE_COMPARISON_SUMMARY, case_comparison_header())
    row = rows[findfirst(row -> row["case_id"] == case_id &&
                                row["conduct"] == CONDUCTIVITY_VARIANT,
                         rows)]

    println("CICE comparison status")
    println("======================")
    println("Current case: $(row["case_id"]) conduct=$(row["conduct"])")
    println("History maps into a ClimaSeaIce state with status $(row["state_replay_status"]).")
    println("Thermodynamic transition status: $(row["case_transition_status"])")
    println("CICE thickness-rate budget: $(row["thickness_rate_budget_status"]) " *
            "(residual $(row["thickness_rate_budget_residual_m"]) m)")
    println("Prognostic ClimaSeaIce thickness replay: $(row["prognostic_thickness_replay_status"]) " *
            "(acceptance ratio $(row["prognostic_thickness_forced_replay_acceptance_ratio"]))")
    println("Source-level prognostic thickness replay: $(row["icepack_matrix_prognostic_thickness_forced_status"]) " *
            "(acceptance ratio $(row["icepack_matrix_prognostic_thickness_forced_replay_acceptance_ratio"]))")
    println("Strict BL99 validation status: $(row["strict_bl99_validation_status"])")
    println("Reset fixed-grid forced replay: $(row["reset_forced_step_status"]) " *
            "(max relative temperature error $(row["reset_forced_step_max_relative_temperature_error"]))")
    println("Reset fixed-grid forced temperature gate: $(row["reset_forced_step_temperature_acceptance_status"]) " *
            "(acceptance ratio $(row["reset_forced_step_temperature_acceptance_ratio"]))")
    println("Icepack-matrix remapped reset replay: $(row["icepack_temperature_matrix_thickness_remapped_reset_status"]) " *
            "(max relative temperature error $(row["icepack_temperature_matrix_thickness_remapped_reset_max_relative_temperature_error"]))")
    println("Icepack-matrix remapped reset temperature gate: $(row["icepack_temperature_matrix_thickness_remapped_reset_temperature_acceptance_status"]) " *
            "(acceptance ratio $(row["icepack_temperature_matrix_thickness_remapped_reset_temperature_acceptance_ratio"]))")
    println("Thickness-remapped reset replay: $(row["thickness_remapped_reset_status"]) " *
            "(max relative temperature error $(row["thickness_remapped_reset_max_relative_temperature_error"]))")
    println("Thickness-remapped reset temperature gate: $(row["thickness_remapped_reset_temperature_acceptance_status"]) " *
            "(acceptance ratio $(row["thickness_remapped_reset_temperature_acceptance_ratio"]))")
    println("CICE-thickness-forced replay: $(row["cice_thickness_forced_status"]) " *
            "(max relative temperature error $(row["cice_thickness_forced_max_relative_temperature_error"]))")
    println("CICE-thickness-forced temperature gate: $(row["cice_thickness_forced_temperature_acceptance_status"]) " *
            "(acceptance ratio $(row["cice_thickness_forced_temperature_acceptance_ratio"]))")
    println("Moving-metric CICE-thickness-forced replay: $(row["moving_metric_cice_thickness_forced_status"]) " *
            "(max relative temperature error $(row["moving_metric_cice_thickness_forced_max_relative_temperature_error"]))")
    println("Moving-metric CICE-thickness-forced temperature gate: $(row["moving_metric_cice_thickness_forced_temperature_acceptance_status"]) " *
            "(acceptance ratio $(row["moving_metric_cice_thickness_forced_temperature_acceptance_ratio"]))")
    println("Icepack-matrix CICE-thickness-forced replay: $(row["icepack_temperature_matrix_cice_thickness_forced_status"]) " *
            "(max relative temperature error $(row["icepack_temperature_matrix_cice_thickness_forced_max_relative_temperature_error"]))")
    println("Icepack-matrix CICE-thickness-forced temperature gate: $(row["icepack_temperature_matrix_cice_thickness_forced_temperature_acceptance_status"]) " *
            "(acceptance ratio $(row["icepack_temperature_matrix_cice_thickness_forced_temperature_acceptance_ratio"]))")
    println("Top-ablation remapped reset replay: $(row["top_ablation_remapped_reset_status"]) " *
            "(max relative temperature error $(row["top_ablation_remapped_reset_max_relative_temperature_error"]))")
    println("Top-ablation remapped reset temperature gate: $(row["top_ablation_remapped_reset_temperature_acceptance_status"]) " *
            "(acceptance ratio $(row["top_ablation_remapped_reset_temperature_acceptance_ratio"]))")
    println("Free-running fixed-grid replay: $(row["free_running_status"]) " *
            "(max relative temperature error $(row["free_running_max_relative_temperature_error"]))")
    println("Free-running fixed-grid temperature gate: $(row["free_running_temperature_acceptance_status"]) " *
            "(acceptance ratio $(row["free_running_temperature_acceptance_ratio"]))")
    println("Aggregate case evidence: $(CASE_COMPARISON_SUMMARY)")
end

function parse_smoke_history!(rows)
    history_dir = get(ENV, "CICE_HISTORY_DIR", DEFAULT_CICE_HISTORY_DIR)

    if !isdir(history_dir)
        println()
        println("CICE history parser")
        println("===================")
        println("History directory not found: $history_dir")
        return rows
    end

    require_ncdump()
    records, _, _ = read_cice_history(history_dir)
    write_cice_summary(CICE_SUMMARY, records)
    print_cice_summary(records, history_dir, CICE_SUMMARY)

    replay_metrics = run_climaseaice_state_replay(records)
    write_replay_summary(CLIMASEAICE_REPLAY_SUMMARY, replay_metrics)
    print_replay_summary(replay_metrics, CLIMASEAICE_REPLAY_SUMMARY)

    forced_replay_metrics = run_fixed_grid_forced_replay(records)
    write_replay_summary(FORCED_REPLAY_SUMMARY, forced_replay_metrics)
    print_forced_replay_summary(forced_replay_metrics, FORCED_REPLAY_SUMMARY)
    upsert_case_comparison_summary!(CASE_COMPARISON_SUMMARY,
                                    history_dir,
                                    records,
                                    replay_metrics,
                                    forced_replay_metrics)

    upsert_dashboard_row!(rows, Dict("gate" => "cice_history_parser",
                                     "status" => "pass",
                                     "evidence" => CICE_SUMMARY))
    upsert_dashboard_row!(rows, Dict("gate" => "climaseaice_state_replay",
                                     "status" => all(metric.status == "pass" for metric in replay_metrics) ? "pass" : "fail",
                                     "evidence" => CLIMASEAICE_REPLAY_SUMMARY))
    upsert_dashboard_row!(rows, Dict("gate" => "fixed_grid_forced_replay",
                                     "status" => "diagnostic",
                                     "evidence" => FORCED_REPLAY_SUMMARY))
    upsert_dashboard_row!(rows, Dict("gate" => "strict_bl99_validation",
                                     "status" => strict_case_comparison_status(CASE_COMPARISON_SUMMARY),
                                     "evidence" => CASE_COMPARISON_SUMMARY))
    for row in rows
        if row["gate"] == "cice_case_comparison"
            row["status"] = "diagnostic"
            row["evidence"] = CASE_COMPARISON_SUMMARY
        end
    end

    write_dashboard(DASHBOARD, rows)
    return rows
end

function main()
    rows = read_dashboard(DASHBOARD)
    history_dir = get(ENV, "CICE_HISTORY_DIR", DEFAULT_CICE_HISTORY_DIR)
    rows = parse_smoke_history!(rows)

    println()
    print_dashboard(rows)

    println()
    if isdir(history_dir)
        print_current_case_status(history_dir)
    else
        println("CICE comparison status")
        println("======================")
        println("CICE history directory not found: $history_dir")
    end
end

main()
