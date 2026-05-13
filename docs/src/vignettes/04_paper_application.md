# Paper Application: Stay-at-Home Orders

```@meta
CurrentModule = MSC
```

Shen, Song, and Abadie apply MSC to estimate the effect of COVID-19
stay-at-home orders on county unemployment in April 2020. Their setup is:

- outcome: county monthly unemployment rate;
- controls: counties in Arkansas, Iowa, Nebraska, North Dakota, South Dakota, Utah, and Wyoming;
- treated units: counties in states with stay-at-home orders;
- post-treatment period: April 2020;
- estimand: county-level ATT for treated counties.

The draft reports 438 control counties, 2,674 treated counties, and an ATT of
about 5.06 percentage points. The text says the data start in January 2010, but
also reports `T0 = 147`; with April 2020 as the post period, `T0 = 147`
corresponds to January 2008 through March 2020. The example script therefore
defaults to January 2008 so that the pre-period count matches the draft.
The replication script also excludes Alaska and the District of Columbia, which
aligns the BLS county file with the paper's reported treated-unit count.

## Data

The public data source is the BLS Local Area Unemployment Statistics flat-file
archive:

```text
https://download.bls.gov/pub/time.series/la/
```

BLS sometimes blocks automated bulk downloads. For that reason, the package
does not download hundreds of megabytes during documentation builds. Instead,
download the county LAUS data file manually from the BLS flat-file archive and
save it as `data/la.data.64.County` from the directory where you run the
example. You can also set `BLS_LAUS_COUNTY_FILE` to an explicit file path.

## Run The Replication Script

From the package root:

```bash
julia --project=. examples/03_covid_sah_orders.jl
```

The script:

1. reads the BLS county unemployment-rate series;
2. keeps complete county histories through April 2020;
3. excludes Alaska and the District of Columbia to match the reported sample;
4. uses the seven no-order states as controls;
5. builds the common-adoption panel;
6. estimates MSC with the application value `lambda = 0.03`;
7. reports ATT, selected `lambda`, and pre-treatment RMSE.

The core estimation call is:

```julia
est = msc_estimate(panel; lambdas = [0.03], standardize = false)
```

Set `MSC_CV=true` to run a cross-validated penalty path instead of the fixed
application value. The example defaults to `MSC_MAX_ITER=100` and
`MSC_TOL=1e-4` so the large county-level application returns in a reasonable
time on a laptop; increase those environment variables for a tighter solver
tolerance.

With the current BLS flat file downloaded under `data/`, the script constructs
the paper-sized panel:

```text
Controls: 438 counties in Arkansas, Iowa, Nebraska, North Dakota, South Dakota, Utah, Wyoming
Treated:  2674 counties
T0:       147 pre-treatment months
ATT = 4.9552
```

That estimate uses `lambda = 0.03` and the default bounded solver settings. The
draft benchmark is about 5.06 percentage points.

## Interpreting Differences

Exact equality with the draft estimate is not guaranteed unless the same BLS
vintage, county inclusion rules, and tuning details are used. LAUS data are
revised, and the paper draft does not include a replication archive. The script
is therefore best read as a transparent replication workflow for the paper's
application, not as a frozen reproduction artifact.
