using MSC
using DataFrames

sim = simulate_msc(ncontrols=12, ntreated=3, T0=40, T1=5, tau=2.0, seed=4)

rows = NamedTuple[]
for i in eachindex(sim.units), t in eachindex(sim.times)
    push!(rows, (unit=sim.units[i], time=sim.times[t], y=sim.Y[i, t], treated=sim.W[i, t]))
end

df = DataFrame(rows)
est = msc_estimate(df, :unit, :time, :y, :treated; lambdas=[0.002], max_iter=1000)

println(est)
println(first(weights_table(est; atol=1e-8), 10))
