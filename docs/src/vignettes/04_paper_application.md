# Paper Application: County Unemployment in April 2020

```@meta
CurrentModule = MSC
```

This vignette is the MSC.jl version of the empirical application in Shen,
Song, and Abadie: the effect of COVID-19 stay-at-home orders on county
unemployment in April 2020.

## Target From The Paper

The paper's application is a high-dimensional synthetic-control problem:

- outcome: county monthly unemployment rate
- controls: counties in Arkansas, Iowa, Nebraska, North Dakota, South Dakota, Utah, and Wyoming
- treated units: counties in states with stay-at-home orders
- post-treatment period: April 2020
- target estimand: ATT across treated counties

The draft reports:

| quantity | paper |
| --- | ---: |
| control counties | 438 |
| treated counties | 2,674 |
| pre-treatment months | 147 |
| April 2020 ATT | about 5.06 percentage points |

The text says the data start in January 2010, but also reports `T0 = 147`.
With April 2020 as the post period, `T0 = 147` corresponds to January 2008
through March 2020. The example script therefore defaults to January 2008.
It also excludes Alaska and the District of Columbia, which aligns the BLS
county file with the paper's reported treated-unit count.

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

## Run The Example

From the package root:

```bash
julia --project=. examples/03_covid_sah_orders.jl
```

The script is [`examples/03_covid_sah_orders.jl`](https://github.com/xiangao/MSC.jl/blob/main/examples/03_covid_sah_orders.jl).
It does the full data construction:

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
application value. The full application has 438 control counties and 2,674
treated counties, so the example defaults to `MSC_MAX_ITER=100` and
`MSC_TOL=1e-4` to keep the script runnable on a laptop.

## What We Replicated

Using the BLS flat file downloaded under `data/`, the example constructs the
same-sized county panel as the paper:

```text
Controls: 438 counties in Arkansas, Iowa, Nebraska, North Dakota, South Dakota, Utah, Wyoming
Treated:  2674 counties
T0:       147 pre-treatment months
```

The script now prints the paper comparison directly:

```text
Paper application replication check
--------------------------------------------------------------------
quantity                             MSC.jl          paper
--------------------------------------------------------------------
control counties                        438            438
treated counties                       2674          2,674
pre-treatment months                    147            147
April 2020 ATT, pp                   4.9552           5.06
--------------------------------------------------------------------
Sample counts match the paper. The ATT is an application replication
using the current BLS flat-file vintage and the script's solver settings.
```

The important distinction is that MSC.jl currently replicates the application
workflow and the reported sample counts. The numerical ATT is close under the
default bounded run, but it is not a frozen exact reproduction of the published
number because the paper does not provide a replication archive with the exact
BLS vintage, county inclusion file, and optimizer/tuning settings.

## Interpreting Differences

If exact equality with 5.06 is the goal, the next step is not more prose in the
vignette. It is to pin the original paper's replication archive or tune the
optimizer against an explicitly stated target. Without that archive, this
example should be read as a transparent application replication: same design,
same sample size, public BLS data, and an openly reported estimate.
