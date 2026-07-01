# MSC.jl — project notes for Claude

Implements the multivariate synthetic control (MSC) estimator of Shen, Song,
and Abadie (2025), "Efficiently Learning Synthetic Control Models for
High-dimensional Disaggregated Data." Fits `Y_pre = X_pre * Theta + E` for
many treated units at once by minimizing a sparse multivariate square-root
Lasso objective:

```
min_Theta (1 / sqrt(T0)) * ||Y_pre - X_pre * Theta||_* + lambda * ||Theta||_1
```

`||.||_*` is the nuclear norm (sum of singular values) of the pre-treatment
residual matrix, so the loss couples all treated units jointly (as opposed to
fitting a separate Lasso per treated unit). `Theta` is `n_controls x
n_treated` and can carry negative weights (unlike classical synthetic
control).

## Main API

- `fit_msc(Xpre, Ypre; lambdas, nlambda, nfolds, ...)` → `MSCFit` — core solver.
- `predict_counterfactual(fit, Xpost)`, `att(fit, Xpost, Ytreated_post)`.
- `msc_estimate(Y, N0, T0; ...)` / `msc_estimate(panel::MSCPanel; ...)` /
  `msc_estimate(table, unit, time, outcome, treatment; ...)` → `MSCEstimate`,
  the matrix / `MSCPanel` / DataFrame entry points.
- `panel_matrices(table, unit, time, outcome, treatment)` — long DataFrame →
  `MSCPanel` (validates common-adoption design: controls never treated,
  treated units switch on together at `T0+1`).
- `simulate_msc(; ncontrols, ntreated, T0, T1, rank, sparsity, tau, noise,
  seed)` — synthetic DGP used throughout tests/vignettes.
- Diagnostics: `effect_matrix`, `effect_curve`, `unit_effects`,
  `weights_table`, `cv_table`, `placebo`, `placebo_distribution`,
  `placebo_se`. Plot recipes via `RecipesBase` for `MSCEstimate` (effect
  curve) and `MSCFit` (objective path or CV-error curve).

## What's where

- `src/MSC.jl` — everything: struct defs, the FISTA-style accelerated
  proximal-gradient solver (`_msrl_single`), lambda-path warm starts
  (`_fit_path`, `_lambda_path`, `_lambda_max`), k-fold CV (`_cv_errors`),
  panel/DataFrame plumbing, and placebo utilities. Single file, ~740 lines.
- `test/runtests.jl` — 4 testsets, 26 `@test`s total.
- `examples/03_covid_sah_orders.jl` — reproduces the paper's county-level
  COVID-SAH-orders application against a BLS LAUS county flat file (not
  bundled; must be downloaded to `data/la.data.64.County`). Sample counts
  match the paper (438 controls, 2,674 treated, 147 pre-periods); the ATT is
  an application replication on the current BLS vintage, not a frozen
  bit-for-bit reproduction.

## Tests

```bash
cd ~/projects/software/MSC.jl && julia --project=. test/runtests.jl
```

26/26 pass, ~20s (verified 2026-07-01).

## Docs

Documenter.jl. `docs/make.jl` builds (`warnonly = true`, `checkdocs = :none`);
`.github/workflows/docs.yml` deploys to `gh-pages` on push to `main`. Pages
include 4 Quarto-authored vignettes under `docs/src/vignettes/` (getting
started, DataFrame panels, diagnostics/placebos, paper application) plus
`reference.md`. Live at <https://xiangao.github.io/MSC.jl/>. Build locally
with `julia --project=docs docs/make.jl`.

## Known warning: "MSC did not converge for lambda=..."

`_msrl_single` runs FISTA (accelerated proximal gradient with backtracking
line search) for up to `max_iter` iterations per lambda, declaring
convergence when the relative change in `Theta` drops below `tol`. `fit_msc`
warns if the *selected* lambda's fit (the CV-chosen one, or the single
lambda passed via `lambdas=`) didn't hit that criterion.

This fires reproducibly in the test suite (all 4 testsets pass regardless)
whenever a test forces `lambdas=[0.002]` together with a low `max_iter`
(1200/800/500 in the "simulation and DataFrame panel API" testset and the
`placebo`/`placebo_distribution` calls that inherit it) — the cap is chosen
for test speed, not for convergence to `tol=1e-7`. It is expected there, not
a bug.

If a user hits this on real data, it means the *chosen* lambda's solver
simply ran out of iterations before hitting the relative-change tolerance —
it is not evidence of divergence or a bad estimate, just an unconfirmed
convergence certificate. What to try, in order of cost:
1. Raise `max_iter` (default 10,000) — cheapest fix, since `_fit_path` warm-starts
   each lambda from the previous one, so later-path fits are usually cheap.
2. Loosen `tol` (default `1e-6`) if the application doesn't need 6 digits of
   agreement in `Theta`.
3. If it persists even with a generous `max_iter`, it usually means `lambda`
   is small relative to the data scale — small `lambda` shrinks the proximal
   step's soft-threshold, so the L1 penalty does little sparsification and
   the nuclear-norm loss dominates; the iterate map is closer to plain
   accelerated gradient descent on a possibly ill-conditioned `X_pre`
   (check `opnorm(Xpre)` / near-collinear controls). Standardizing
   (`standardize=true`) or using a coarser/higher `lambda_min_ratio` lambda
   grid (so the smallest lambda tried isn't as extreme) often resolves it
   faster than just raising `max_iter` further.
