# Examples

```@meta
CurrentModule = MSC
```

Runnable scripts live in the repository's `examples/` directory. These mirror
the package vignettes and are useful when you want a plain Julia script instead
of a documentation page.

## Matrix Workflow

```bash
julia --project=. examples/01_matrix_workflow.jl
```

This script simulates a common-adoption panel, estimates MSC from the panel
matrix, prints the ATT, reports treated-unit effects, and summarizes the
pre-treatment fit.

Source:
[`examples/01_matrix_workflow.jl`](https://github.com/xiangao/MSC.jl/blob/main/examples/01_matrix_workflow.jl)

## DataFrame Workflow

```bash
julia --project=. examples/02_dataframe_workflow.jl
```

This script converts a long `DataFrame` with `unit`, `time`, `y`, and
`treated` columns into the MSC panel format, estimates the ATT, and prints the
largest nonzero control-treated weights.

Source:
[`examples/02_dataframe_workflow.jl`](https://github.com/xiangao/MSC.jl/blob/main/examples/02_dataframe_workflow.jl)

## Paper Application Scaffold

```bash
julia --project=. examples/03_covid_sah_orders.jl
```

This script implements the package's replication scaffold for Shen, Song, and
Abadie's COVID-19 stay-at-home orders application. It expects a local BLS LAUS
county flat file at `data/la.data.64.County`, or an explicit
`BLS_LAUS_COUNTY_FILE`, because BLS bulk-download access is sometimes blocked
for scripted requests.

Source:
[`examples/03_covid_sah_orders.jl`](https://github.com/xiangao/MSC.jl/blob/main/examples/03_covid_sah_orders.jl)
