# MSC.jl

`MSC.jl` implements multivariate synthetic control for high-dimensional
disaggregated panels with many treated units. It is inspired by Shen, Song, and
Abadie, "Efficiently Learning Synthetic Control Models for High-dimensional
Disaggregated Data" (2025).

The estimator fits

```text
Y_pre = X_pre * Theta + E
```

where `X_pre` is a `T0 x n_controls` matrix and `Y_pre` is a
`T0 x n_treated` matrix. The sparse multivariate square-root Lasso objective is

```text
min_Theta (1 / sqrt(T0)) * ||Y_pre - X_pre * Theta||_* + lambda * ||Theta||_1
```

## Installation

```julia
using Pkg
Pkg.develop(path = "/home/xao/projects/software/MSC.jl")
```

## Matrix API

```julia
using MSC

fit = fit_msc(Xpre, Ypre; nlambda = 50, nfolds = 5)
y0 = predict_counterfactual(fit, Xpost)
tau = att(fit, Xpost, Ytreated_post)
```

For a panel matrix with controls in the first `N0` rows and treated units in the
remaining rows:

```julia
est = msc_estimate(Y, N0, T0; nlambda = 50, nfolds = 5)
est.estimate
```

## DataFrame API

```julia
est = msc_estimate(df, :unit, :time, :outcome, :treated)
weights_table(est)
effect_curve(est)
```

## Vignettes and Examples

Documentation vignettes:

- [Getting Started](docs/src/vignettes/01_getting_started.md)
- [DataFrame Panels](docs/src/vignettes/02_dataframe_panels.md)
- [Diagnostics and Placebos](docs/src/vignettes/03_diagnostics_placebos.md)
- [Paper Application](docs/src/vignettes/04_paper_application.md)

Runnable examples:

- [Matrix workflow](examples/01_matrix_workflow.jl)
- [DataFrame workflow](examples/02_dataframe_workflow.jl)
- [COVID stay-at-home orders application](examples/03_covid_sah_orders.jl)

## Development

```julia
using Pkg
Pkg.test("MSC")
```

Build docs locally:

```bash
julia --project=docs docs/make.jl
```

Run the paper-application scaffold with a local BLS LAUS county flat file:

```bash
BLS_LAUS_COUNTY_FILE=/path/to/la.data.64.County \
    julia --project=. examples/03_covid_sah_orders.jl
```
