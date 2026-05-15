using MSC
using DataFrames
using LinearAlgebra
using Random
using Statistics
using Test

@testset "MSC core" begin
    rng = MersenneTwister(42)
    T0 = 50
    T1 = 4
    ncontrols = 12
    ntreated = 3

    Xpre = randn(rng, T0, ncontrols)
    theta = zeros(ncontrols, ntreated)
    theta[1, 1] = 0.7
    theta[2, 1] = 0.2
    theta[3, 2] = 0.5
    theta[4, 2] = -0.3
    theta[5, 3] = 1.0
    Ypre = Xpre * theta .+ 0.02 .* randn(rng, T0, ntreated)

    fit = fit_msc(Xpre, Ypre; lambdas=[0.005], max_iter=2000, tol=1e-7)
    @test size(coef(fit)) == (ncontrols, ntreated)
    @test objective(Xpre, Ypre, coef(fit), fit.lambda) < objective(Xpre, Ypre, zeros(ncontrols, ntreated), fit.lambda)

    Xpost = randn(rng, T1, ncontrols)
    y0 = predict_counterfactual(fit, Xpost)
    @test size(y0) == (T1, ntreated)

    effect = 2.0
    Ypost = Xpost * theta .+ effect
    @test isapprox(att(fit, Xpost, Ypost), effect; atol=0.15)
end

@testset "panel wrapper and CV" begin
    rng = MersenneTwister(7)
    T0 = 30
    T1 = 2
    ncontrols = 8
    ntreated = 2

    X = randn(rng, T0 + T1, ncontrols)
    theta = zeros(ncontrols, ntreated)
    theta[1, :] .= [0.4, -0.2]
    theta[2, :] .= [0.3, 0.5]
    Y0 = X * theta
    tau = 1.25
    Ytreated = copy(Y0)
    Ytreated[(T0 + 1):end, :] .+= tau
    panel = vcat(X', Ytreated')

    est = msc_estimate(panel, ncontrols, T0; nlambda=5, nfolds=3, rng=rng, max_iter=1000, tol=1e-6)
    @test est isa MSCEstimate
    @test length(est.fit.cv_errors) == 5
    @test isapprox(est.estimate, tau; atol=0.2)
    @test size(effect_matrix(est)) == (T1, ntreated)
    @test length(effect_curve(est)) == T1
    @test length(unit_effects(est)) == ntreated
    @test nrow(weights_table(est; atol=1e-8)) > 0
    @test nrow(cv_table(est.fit)) == 5
end

@testset "convergence, cv_table, and intercept" begin
    rng = MersenneTwister(17)
    T0 = 40
    ncontrols = 6
    ntreated = 2
    Xpre = randn(rng, T0, ncontrols)
    theta_true = zeros(ncontrols, ntreated)
    theta_true[1, 1] = 0.6; theta_true[2, 2] = 0.8
    Ypre = Xpre * theta_true .+ 0.01 .* randn(rng, T0, ntreated)

    # Convergence flag should be true for an easy, low-noise problem
    fit = fit_msc(Xpre, Ypre; lambdas=[0.003], max_iter=5000, tol=1e-7)
    @test fit.converged

    # cv_table: selected column should have exactly one true, at the correct index
    fit_cv = fit_msc(Xpre, Ypre; nlambda=8, nfolds=4, rng=rng, max_iter=2000, tol=1e-6)
    tbl = cv_table(fit_cv)
    @test sum(tbl.selected) == 1
    @test tbl.selected[fit_cv.lambda_idx]
    @test tbl.lambda[fit_cv.lambda_idx] == fit_cv.lambda

    # Intercept round-trip: predict_counterfactual on training data should match
    # the mean of Ypre when Xpre is the design matrix (i.e., intercept is applied correctly)
    yhat_pre = predict_counterfactual(fit, Xpre)
    @test size(yhat_pre) == (T0, ntreated)
    # Residual mean should be small (intercept centers the prediction)
    @test all(abs.(vec(mean(Ypre .- yhat_pre; dims=1))) .< 0.05)
end

@testset "simulation and DataFrame panel API" begin
    sim = simulate_msc(ncontrols=10, ntreated=3, T0=35, T1=4, tau=1.5, seed=11, noise=0.01)
    est = msc_estimate(sim.Y, sim.N0, sim.T0; lambdas=[0.002], max_iter=1200, tol=1e-7)
    @test isapprox(est.estimate, 1.5; atol=0.2)

    units = sim.units
    times = sim.times
    rows = NamedTuple[]
    for i in eachindex(units), t in eachindex(times)
        push!(
            rows,
            (
                unit=units[i],
                time=times[t],
                y=sim.Y[i, t],
                treated=sim.W[i, t],
            ),
        )
    end
    df = DataFrame(rows)
    panel = panel_matrices(df, :unit, :time, :y, :treated)
    @test panel isa MSCPanel
    @test panel.N0 == sim.N0
    @test panel.T0 == sim.T0

    est_df = msc_estimate(df, :unit, :time, :y, :treated; lambdas=[0.002], max_iter=1200, tol=1e-7)
    @test isapprox(est_df.estimate, est.estimate; atol=1e-8)

    plac = placebo(est_df; T0_placebo=20, lambdas=[0.002], max_iter=800, tol=1e-6)
    @test plac isa MSCEstimate
    dist = placebo_distribution(est_df; replications=3, ntreated=2, lambdas=[0.002], max_iter=500, tol=1e-6, rng=MersenneTwister(9))
    @test length(dist) == 3
    @test isfinite(placebo_se(est_df; replications=3, ntreated=2, lambdas=[0.002], max_iter=500, tol=1e-6, rng=MersenneTwister(9)))
end
