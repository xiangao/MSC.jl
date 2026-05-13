# MSC.jl

`MSC.jl` implements multivariate synthetic control for high-dimensional
disaggregated panels. It is designed for settings with many treated units and many
controls where fitting one synthetic-control model per treated unit is expensive.

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

## Vignettes

- [Getting Started](vignettes/01_getting_started.md): matrix workflow with simulated data.
- [DataFrame Panels](vignettes/02_dataframe_panels.md): long-panel conversion and direct table estimation.
- [Diagnostics and Placebos](vignettes/03_diagnostics_placebos.md): pre-fit checks, CV output, and placebo routines.
- [Paper Application](vignettes/04_paper_application.md): replication scaffold for the COVID stay-at-home order example.

## Examples

Runnable scripts live in the repository:

- [`examples/01_matrix_workflow.jl`](https://github.com/xiangao/MSC.jl/blob/main/examples/01_matrix_workflow.jl)
- [`examples/02_dataframe_workflow.jl`](https://github.com/xiangao/MSC.jl/blob/main/examples/02_dataframe_workflow.jl)
- [`examples/03_covid_sah_orders.jl`](https://github.com/xiangao/MSC.jl/blob/main/examples/03_covid_sah_orders.jl)
