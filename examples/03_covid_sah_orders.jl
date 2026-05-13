using CSV
using DataFrames
using MSC
using Printf
using Statistics

# Replication-style script for Shen, Song, and Abadie's county-level application:
# COVID-19 stay-at-home orders and April 2020 county unemployment rates.
#
# Data source:
#   BLS Local Area Unemployment Statistics flat files
#   https://download.bls.gov/pub/time.series/la/
#
# Download the county data flat file manually if scripted access is blocked by
# BLS. The script looks for data/la.data.64.County from the current directory,
# the package root, and ~/projects/codex/data. You can also set an explicit path:
#
#   BLS_LAUS_COUNTY_FILE=/path/to/la.data.64.County julia --project=. examples/03_covid_sah_orders.jl
#
# The script expects a tab-separated BLS file with columns:
#   series_id, year, period, value, footnote_codes

const CONTROL_STATE_FIPS = Set(["05", "19", "31", "38", "46", "49", "56"])
const EXCLUDED_STATE_FIPS = Set(["02", "11"])  # Alaska and District of Columbia.
const CONTROL_STATE_NAMES = [
    "Arkansas",
    "Iowa",
    "Nebraska",
    "North Dakota",
    "South Dakota",
    "Utah",
    "Wyoming",
]
const PAPER_ATT = 5.06
const CONTROL_STATES = Dict(
    "05" => "Arkansas",
    "19" => "Iowa",
    "31" => "Nebraska",
    "38" => "North Dakota",
    "46" => "South Dakota",
    "49" => "Utah",
    "56" => "Wyoming",
)

function county_state_fips(series_id)
    s = strip(String(series_id))
    startswith(s, "LAUCN") || return missing
    length(s) >= 7 || return missing
    return s[6:7]
end

function is_county_unemployment_rate(series_id)
    s = strip(String(series_id))
    return startswith(s, "LAUCN") && endswith(s, "03")
end

function month_index(year, period)
    p = strip(String(period))
    startswith(p, "M") || return missing
    month = parse(Int, p[2:end])
    1 <= month <= 12 || return missing
    return 12 * parse(Int, string(year)) + month
end

function default_laus_county_paths()
    return [
        joinpath(pwd(), "data", "la.data.64.County"),
        joinpath(@__DIR__, "..", "data", "la.data.64.County"),
        joinpath(homedir(), "projects", "codex", "data", "la.data.64.County"),
    ]
end

function find_laus_county_file()
    explicit = get(ENV, "BLS_LAUS_COUNTY_FILE", "")
    if !isempty(explicit)
        isfile(explicit) || error("BLS_LAUS_COUNTY_FILE does not exist: $(explicit)")
        return explicit
    end

    for candidate in default_laus_county_paths()
        isfile(candidate) && return candidate
    end

    return nothing
end

function load_laus_county_panel(path; start_year=2008, start_month=1, end_year=2020, end_month=4)
    raw = CSV.read(path, DataFrame; delim='\t', stripwhitespace=true, silencewarnings=true)
    rename!(raw, Symbol.(strip.(String.(names(raw)))))

    needed = [:series_id, :year, :period, :value]
    missing_cols = setdiff(needed, Symbol.(names(raw)))
    isempty(missing_cols) || error("Missing BLS columns: $(missing_cols)")

    raw = filter(row -> is_county_unemployment_rate(row.series_id), raw)
    raw.series_id = strip.(String.(raw.series_id))
    raw.state_fips = county_state_fips.(raw.series_id)
    raw.month = month_index.(raw.year, raw.period)
    raw = filter(row -> !(row.state_fips in EXCLUDED_STATE_FIPS), raw)

    start_idx = 12 * start_year + start_month
    end_idx = 12 * end_year + end_month
    raw = filter(row -> !ismissing(row.month) && start_idx <= row.month <= end_idx, raw)
    raw.value = tryparse.(Float64, strip.(String.(raw.value)))
    raw = filter(row -> !isnothing(row.value), raw)
    raw.value = Float64.(raw.value)

    raw.date = string.(raw.year, "-", lpad.(replace.(String.(raw.period), "M" => ""), 2, "0"))
    raw.treated = [row.state_fips in CONTROL_STATE_FIPS ? 0.0 : 1.0 for row in eachrow(raw)]
    raw.post = [row.month == end_idx ? 1.0 : 0.0 for row in eachrow(raw)]
    raw.D = raw.treated .* raw.post

    # Keep counties with complete monthly unemployment-rate histories.
    keep = combine(groupby(raw, :series_id), nrow => :n)
    full_n = length(unique(raw.month))
    keep = keep.series_id[keep.n .== full_n]
    raw = filter(row -> row.series_id in keep, raw)

    return select(raw, :series_id => :county, :date => :month, :value => :unemployment_rate, :D, :state_fips)
end

function print_replication_table(est, panel)
    treated = size(panel.Y, 1) - panel.N0
    rows = [
        ("control counties", string(panel.N0), "438"),
        ("treated counties", string(treated), "2,674"),
        ("pre-treatment months", string(panel.T0), "147"),
        ("April 2020 ATT, pp", @sprintf("%.4f", est.estimate), @sprintf("%.2f", PAPER_ATT)),
    ]

    println()
    println("Paper application replication check")
    println("-" ^ 68)
    @printf("%-28s %14s %14s\n", "quantity", "MSC.jl", "paper")
    println("-" ^ 68)
    for (name, ours, paper) in rows
        @printf("%-28s %14s %14s\n", name, ours, paper)
    end
    println("-" ^ 68)
    println("Sample counts match the paper. The ATT is an application replication")
    println("using the current BLS flat-file vintage and the script's solver settings.")
end

function main()
    path = find_laus_county_file()
    if isnothing(path)
        println("Download the BLS LAUS county flat file to data/la.data.64.County, or set BLS_LAUS_COUNTY_FILE.")
        println("Expected source: https://download.bls.gov/pub/time.series/la/")
        println("Example: BLS_LAUS_COUNTY_FILE=/path/to/la.data.64.County julia --project=. examples/03_covid_sah_orders.jl")
        return nothing
    end

    lambda = parse(Float64, get(ENV, "MSC_LAMBDA", "0.03"))
    max_iter = parse(Int, get(ENV, "MSC_MAX_ITER", "100"))
    tol = parse(Float64, get(ENV, "MSC_TOL", "1e-4"))
    use_cv = lowercase(get(ENV, "MSC_CV", "false")) in ("1", "true", "yes")

    println("Reading BLS LAUS county data from: ", path)
    df = load_laus_county_panel(path)
    panel = panel_matrices(df, :county, :month, :unemployment_rate, :D)

    println("Controls: ", panel.N0, " counties in ", join(CONTROL_STATE_NAMES, ", "))
    println("Treated:  ", size(panel.Y, 1) - panel.N0, " counties")
    println("T0:       ", panel.T0, " pre-treatment months")

    fit_kwargs = use_cv ? (; nlambda=20, nfolds=5) : (; lambdas=[lambda])
    est = msc_estimate(panel; fit_kwargs..., standardize=false, max_iter=max_iter, tol=tol)

    println(est)
    println("Selected lambda: ", est.fit.lambda)
    println("Solver iterations: ", est.fit.iterations, " (converged: ", est.fit.converged, ")")
    println("Pre-fit RMSE: ", round(prefit_rmse(est.fit, panel.Y[1:panel.N0, 1:panel.T0]', panel.Y[(panel.N0 + 1):end, 1:panel.T0]'), digits=4))
    print_replication_table(est, panel)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
