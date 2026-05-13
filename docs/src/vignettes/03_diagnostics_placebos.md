# Diagnostics and Placebos

```@meta
CurrentModule = MSC
```

MSC is useful only when the controls can reconstruct treated-unit pre-treatment
paths. Start with pre-fit diagnostics.

```@example diagnostics
using MSC
using Statistics

sim = simulate_msc(ncontrols = 25, ntreated = 5, T0 = 55, T1 = 6, tau = 1.0, seed = 5)
est = msc_estimate(sim.Y, sim.N0, sim.T0; nlambda = 8, nfolds = 4, max_iter = 900)

Xpre = sim.Y[1:sim.N0, 1:sim.T0]'
Ypre = sim.Y[(sim.N0 + 1):end, 1:sim.T0]'

(
    att = round(est.estimate, digits = 3),
    pre_rmse = round(prefit_rmse(est.fit, Xpre, Ypre), digits = 3),
    selected_lambda = est.fit.lambda,
)
```

Cross-validation output is stored as a table:

```@example diagnostics
cv_table(est.fit)
```

One time-placebo moves the intervention date into the original pre-treatment
period:

```@example diagnostics
pl = placebo(est; T0_placebo = 30, lambdas = [est.fit.lambda], max_iter = 700)
round(pl.estimate, digits = 3)
```

A placebo distribution randomly assigns control units to placebo treatment groups:

```@example diagnostics
dist = placebo_distribution(
    est;
    replications = 10,
    ntreated = 3,
    lambdas = [est.fit.lambda],
    max_iter = 500,
)

(mean = round(mean(dist), digits = 3), se = round(std(dist), digits = 3))
```
