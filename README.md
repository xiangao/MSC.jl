# MSC.jl

`MSC.jl` implements multivariate synthetic control for panels with many treated
units. The estimator follows Shen, Song, and Abadie, "Efficiently Learning
Synthetic Control Models for High-dimensional Disaggregated Data" (2025).

Documentation, vignettes, and examples: https://xiangao.github.io/MSC.jl/

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

## Paper Application

The repository includes a runnable version of the county-unemployment
application from Shen, Song, and Abadie:

```bash
julia --project=. examples/03_covid_sah_orders.jl
```

After downloading the BLS LAUS county flat file to `data/la.data.64.County`,
the script constructs the paper-sized panel and prints the comparison:

```text
quantity                             MSC.jl          paper
control counties                        438            438
treated counties                       2674          2,674
pre-treatment months                    147            147
April 2020 ATT, pp                   4.9552           5.06
```

The sample counts match the paper. The ATT is an application replication using
the current BLS flat-file vintage and the script's solver settings, not a frozen
bit-for-bit reproduction of the paper's private run.

## Development

```julia
using Pkg
Pkg.test("MSC")
```

Build docs locally:

```bash
julia --project=docs docs/make.jl
```

Run the paper-application scaffold after downloading the BLS LAUS county flat
file to `data/la.data.64.County`, or set an explicit path:

```bash
julia --project=. examples/03_covid_sah_orders.jl
```
