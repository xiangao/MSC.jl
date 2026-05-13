using CSV
using DataFrames
using MSC
using Statistics

# Replication-style script for Shen, Song, and Abadie's county-level application:
# COVID-19 stay-at-home orders and April 2020 county unemployment rates.
#
# Data source:
#   BLS Local Area Unemployment Statistics flat files
#   https://download.bls.gov/pub/time.series/la/
#
# Download the county data flat file manually if scripted access is blocked by BLS,
# then run:
#
#   BLS_LAUS_COUNTY_FILE=/path/to/la.data.64.County julia --project=. examples/03_covid_sah_orders.jl
#
# The script expects a tab-separated BLS file with columns:
#   series_id, year, period, value, footnote_codes

const CONTROL_STATE_FIPS = Set(["05", "19", "31", "38", "46", "49", "56"])
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
    s = String(series_id)
    startswith(s, "LAUCN") || return missing
    length(s) >= 7 || return missing
    return s[6:7]
end

function is_county_unemployment_rate(series_id)
    s = String(series_id)
    return startswith(s, "LAUCN") && endswith(s, "03")
end

function month_index(year, period)
    p = String(period)
    startswith(p, "M") || return missing
    month = parse(Int, p[2:end])
    1 <= month <= 12 || return missing
    return 12 * parse(Int, string(year)) + month
end

function load_laus_county_panel(path; start_year=2008, start_month=1, end_year=2020, end_month=4)
    raw = CSV.read(path, DataFrame; delim='\t', ignorerepeated=true, stripwhitespace=true)
    rename!(raw, Symbol.(strip.(String.(names(raw)))))

    needed = [:series_id, :year, :period, :value]
    missing_cols = setdiff(needed, Symbol.(names(raw)))
    isempty(missing_cols) || error("Missing BLS columns: $(missing_cols)")

    raw = filter(row -> is_county_unemployment_rate(row.series_id), raw)
    raw.state_fips = county_state_fips.(raw.series_id)
    raw.month = month_index.(raw.year, raw.period)

    start_idx = 12 * start_year + start_month
    end_idx = 12 * end_year + end_month
    raw = filter(row -> !ismissing(row.month) && start_idx <= row.month <= end_idx, raw)
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

function main()
    path = get(ENV, "BLS_LAUS_COUNTY_FILE", "")
    if isempty(path) || !isfile(path)
        println("Set BLS_LAUS_COUNTY_FILE to a local BLS LAUS county flat file.")
        println("Expected source: https://download.bls.gov/pub/time.series/la/")
        println("Example: BLS_LAUS_COUNTY_FILE=/path/to/la.data.64.County julia --project=. examples/03_covid_sah_orders.jl")
        return nothing
    end

    df = load_laus_county_panel(path)
    panel = panel_matrices(df, :county, :month, :unemployment_rate, :D)

    println("Controls: ", panel.N0, " counties in ", join(values(CONTROL_STATES), ", "))
    println("Treated:  ", size(panel.Y, 1) - panel.N0, " counties")
    println("T0:       ", panel.T0, " pre-treatment months")

    est = msc_estimate(
        panel;
        nlambda=20,
        nfolds=5,
        standardize=false,
        max_iter=1500,
        tol=1e-6,
    )

    println(est)
    println("Selected lambda: ", est.fit.lambda)
    println("Pre-fit RMSE: ", round(prefit_rmse(est.fit, panel.Y[1:panel.N0, 1:panel.T0]', panel.Y[(panel.N0 + 1):end, 1:panel.T0]'), digits=4))
    println("Paper benchmark reported in the draft: about 5.06 percentage points.")
end

main()
