# DataFrame Panels

`MSC.jl` includes a small long-panel helper for common-adoption designs. Treatment
must be binary, controls must never be treated, and treated units must remain
treated after the first post-treatment period.

```@example dataframe_panels
using MSC
using DataFrames

sim = simulate_msc(ncontrols = 12, ntreated = 3, T0 = 40, T1 = 5, tau = 2.0, seed = 4)

rows = NamedTuple[]
for i in eachindex(sim.units), t in eachindex(sim.times)
    push!(
        rows,
        (
            unit = sim.units[i],
            time = sim.times[t],
            outcome = sim.Y[i, t],
            treated = sim.W[i, t],
        ),
    )
end

df = DataFrame(rows)
first(df, 5)
```

Convert the long table to matrix form:

```@example dataframe_panels
panel = panel_matrices(df, :unit, :time, :outcome, :treated)

(size(panel.Y), panel.N0, panel.T0)
```

Estimate directly from the table:

```@example dataframe_panels
est = msc_estimate(
    df,
    :unit,
    :time,
    :outcome,
    :treated;
    lambdas = [0.002],
    max_iter = 1000,
    tol = 1e-6,
)

round(est.estimate, digits = 3)
```

The table labels are preserved for diagnostics:

```@example dataframe_panels
first(weights_table(est; atol = 1e-8), 10)
```
