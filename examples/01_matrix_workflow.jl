using MSC
using Statistics

sim = simulate_msc(ncontrols=30, ntreated=6, T0=60, T1=8, tau=1.5, seed=11)

est = msc_estimate(sim.Y, sim.N0, sim.T0; nlambda=10, nfolds=5, max_iter=1000)

println(est)
println("Unit effects: ", round.(unit_effects(est), digits=3))
println("Pre-fit RMSE: ", round(prefit_rmse(est.fit, sim.Y[1:sim.N0, 1:sim.T0]', sim.Y[(sim.N0 + 1):end, 1:sim.T0]'), digits=3))
