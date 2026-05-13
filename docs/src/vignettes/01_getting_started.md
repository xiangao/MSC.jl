# Getting Started

```@meta
CurrentModule = MSC
```

This vignette uses the matrix API. Rows of `Y` are units, columns are time
periods, controls come first, and treated units come last.

```@example getting_started
using MSC
using Statistics

sim = simulate_msc(
    ncontrols = 30,
    ntreated = 6,
    T0 = 60,
    T1 = 8,
    tau = 1.5,
    noise = 0.03,
    seed = 11,
)

est = msc_estimate(
    sim.Y,
    sim.N0,
    sim.T0;
    nlambda = 10,
    nfolds = 5,
    max_iter = 1000,
    tol = 1e-6,
)

round(est.estimate, digits = 3)
```

`effect_matrix` returns post-treatment effects by period and treated unit. The
ATT is the average of this matrix.

```@example getting_started
effects = effect_matrix(est)
(
    size = size(effects),
    att = round(mean(effects), digits = 3),
    by_unit = round.(unit_effects(est), digits = 3),
)
```

The fitted weight matrix is sparse by construction:

```@example getting_started
first(weights_table(est; atol = 1e-8), 8)
```

Plotting recipes are available when `Plots` is loaded:

```@example getting_started
using Plots

plt = plot(
    est;
    title = "MSC average effect",
    xlabel = "Post-treatment period",
    ylabel = "Effect",
)

savefig(plt, "getting-started-msc.svg") # hide
nothing # hide
```

![](getting-started-msc.svg)
