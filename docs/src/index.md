# MSC.jl

`MSC.jl` implements multivariate synthetic control for high-dimensional
disaggregated panels. It is designed for settings with many treated units and many
controls where fitting one synthetic-control model per treated unit is expensive.

The package includes:

- `fit_msc` and `predict_counterfactual` for matrix-level MSC estimation.
- `msc_estimate` for panel-matrix and long-table ATT estimation.
- `panel_matrices` for balanced common-adoption panel conversion.
- Cross-validated regularization over the MSC penalty path.
- Treated-unit counterfactuals, effect curves, unit effects, and weight tables.
- Time-placebo and control-unit placebo routines.
- A lightweight `RecipesBase` plotting recipe.
- Synthetic data generators and a paper-application replication scaffold.

## Where To Start

- Use [Getting Started](vignettes/01_getting_started.md) for the matrix API.
- Use [DataFrame Panels](vignettes/02_dataframe_panels.md) for long panel data.
- Use [Paper Application](vignettes/04_paper_application.md) for the
  county-unemployment application from Shen, Song, and Abadie.
- Use [Examples](examples.md) for runnable scripts.

The paper application script uses the local BLS data file and prints the
replication comparison directly:

```text
control counties       MSC.jl 438     paper 438
treated counties       MSC.jl 2674    paper 2,674
pre-treatment months   MSC.jl 147     paper 147
April 2020 ATT, pp     MSC.jl 4.9552  paper 5.06
```

The estimator fits

```text
Y_pre = X_pre * Theta + E
```

with the sparse multivariate square-root Lasso objective

```text
min_Theta (1 / sqrt(T0)) * ||Y_pre - X_pre * Theta||_* + lambda * ||Theta||_1
```

where `X_pre` is `T0 x n_controls`, `Y_pre` is `T0 x n_treated`, and `Theta`
maps control outcomes to treated-unit counterfactuals.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/xiangao/MSC.jl")
```

For local development:

```julia
using Pkg
Pkg.develop(path = "/home/xao/projects/software/MSC.jl")
```

## Quick Start

```@example quickstart
using MSC
using Statistics

sim = simulate_msc(ncontrols = 20, ntreated = 5, T0 = 50, T1 = 6, tau = 1.2, seed = 1)

est = msc_estimate(
    sim.Y,
    sim.N0,
    sim.T0;
    nlambda = 8,
    nfolds = 4,
    max_iter = 800,
    tol = 1e-6,
)

round(est.estimate, digits = 3)
```

The fitted object stores treated-unit counterfactuals and diagnostics:

```@example quickstart
(
    effects = round.(unit_effects(est), digits = 3),
    prefit = round(prefit_rmse(est.fit, sim.Y[1:sim.N0, 1:sim.T0]', sim.Y[(sim.N0 + 1):end, 1:sim.T0]'), digits = 3),
    nonzero_weights = size(weights_table(est; atol = 1e-8), 1),
)
```

See the vignettes and examples pages for complete workflows.
