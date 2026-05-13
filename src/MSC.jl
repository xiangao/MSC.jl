module MSC

using DataFrames
using LinearAlgebra
using Random
using RecipesBase
using Statistics
using Tables

export MSCFit,
    MSCEstimate,
    MSCPanel,
    fit_msc,
    predict_counterfactual,
    att,
    msc_estimate,
    panel_matrices,
    simulate_msc,
    coef,
    objective,
    fitted_pre,
    residuals,
    prefit_rmse,
    effect_matrix,
    effect_curve,
    unit_effects,
    weights_table,
    cv_table,
    placebo,
    placebo_distribution,
    placebo_se

"""
    MSCFit

Fitted multivariate synthetic-control model. `theta` is an `n_controls x n_treated`
matrix mapping control outcomes to treated counterfactual outcomes.
"""
struct MSCFit
    theta::Matrix{Float64}
    lambda::Float64
    lambdas::Vector{Float64}
    cv_errors::Union{Nothing,Vector{Float64}}
    objective_values::Vector{Float64}
    iterations::Int
    converged::Bool
    x_mean::Vector{Float64}
    x_scale::Vector{Float64}
    y_mean::Vector{Float64}
    standardize::Bool
end

"""
    MSCEstimate

Point estimate and fitted object for an MSC ATT analysis.
Rows of `counterfactual` are post-treatment periods and columns are treated units.
"""
struct MSCEstimate
    estimate::Float64
    fit::MSCFit
    counterfactual::Matrix{Float64}
    treated_post::Matrix{Float64}
    Y::Matrix{Float64}
    N0::Int
    T0::Int
    units::Vector{Any}
    times::Vector{Any}
end

"""
    MSCPanel

Balanced panel in matrix form. Rows `1:N0` are controls, rows `N0+1:end`
are treated, columns `1:T0` are pre-treatment periods, and remaining columns
are post-treatment periods.
"""
struct MSCPanel
    Y::Matrix{Float64}
    N0::Int
    T0::Int
    W::Matrix{Float64}
    units::Vector{Any}
    times::Vector{Any}
end

Base.show(io::IO, fit::MSCFit) = print(
    io,
    "MSCFit(lambda=$(round(fit.lambda, sigdigits=4)), nonzero=$(count(!iszero, fit.theta)), ",
    "iterations=$(fit.iterations), converged=$(fit.converged))",
)

function Base.show(io::IO, est::MSCEstimate)
    println(io, "MSC Estimate")
    println(io, "-" ^ 28)
    println(io, "  ATT = ", round(est.estimate, digits=4))
    println(io, "  lambda = ", round(est.fit.lambda, sigdigits=4))
    println(io, "  nonzero weights = ", count(!iszero, est.fit.theta))
    println(io, "  treated units = ", size(est.treated_post, 2))
    println(io, "  post periods = ", size(est.treated_post, 1))
end

"""
    coef(fit)

Return the `n_controls x n_treated` MSC weight matrix.
"""
coef(fit::MSCFit) = fit.theta

soft_threshold(x::Real, t::Real) = sign(x) * max(abs(x) - t, zero(float(x)))
soft_threshold(A::AbstractMatrix{<:Real}, t::Real) = soft_threshold.(A, t)

function nuclear_norm(A::AbstractMatrix{<:Real})
    isempty(A) && return 0.0
    return sum(svdvals(Matrix{Float64}(A)))
end

"""
    objective(X, Y, theta, lambda)

Evaluate the MSC objective
`(1 / sqrt(size(X, 1))) * ||Y - X * theta||_* + lambda * ||theta||_1`.
"""
function objective(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real},
                   theta::AbstractMatrix{<:Real}, lambda::Real)
    nobs = size(X, 1)
    return nuclear_norm(Y - X * theta) / sqrt(nobs) + lambda * sum(abs, theta)
end

function _default_lambdas(nlambda::Int, lambda_min_ratio::Float64)
    nlambda >= 1 || throw(ArgumentError("nlambda must be positive"))
    0 < lambda_min_ratio <= 1 || throw(ArgumentError("lambda_min_ratio must be in (0, 1]"))
end

function _center_scale(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real}, standardize::Bool)
    Xf = Matrix{Float64}(X)
    Yf = Matrix{Float64}(Y)
    x_mean = vec(mean(Xf; dims=1))
    y_mean = vec(mean(Yf; dims=1))
    Xc = Xf .- reshape(x_mean, 1, :)
    Yc = Yf .- reshape(y_mean, 1, :)
    if standardize
        x_scale = vec(std(Xf; dims=1, corrected=true))
        x_scale .= ifelse.(x_scale .> 0, x_scale, 1.0)
        Xc ./= reshape(x_scale, 1, :)
    else
        x_scale = ones(size(Xf, 2))
    end
    return Xc, Yc, x_mean, x_scale, y_mean
end

function _lambda_max(X::AbstractMatrix{Float64}, Y::AbstractMatrix{Float64})
    nobs = size(X, 1)
    F = svd(Y; full=false)
    isempty(F.S) && return 1.0
    G = X' * (F.U * F.Vt) / sqrt(nobs)
    lm = maximum(abs, G)
    return isfinite(lm) && lm > 0 ? 1.01 * lm : 1.0
end

function _lambda_path(X::AbstractMatrix{Float64}, Y::AbstractMatrix{Float64};
                      nlambda::Int=50, lambda_min_ratio::Float64=1e-2)
    _default_lambdas(nlambda, lambda_min_ratio)
    lambda_max = _lambda_max(X, Y)
    nlambda == 1 && return [lambda_max]
    lambda_min = lambda_min_ratio * lambda_max
    return exp.(range(log(lambda_max), log(lambda_min); length=nlambda))
end

"""
    fitted_pre(fit, Xpre)

Return fitted treated-unit outcomes in the pre-treatment period.
"""
function fitted_pre(fit::MSCFit, Xpre::AbstractMatrix{<:Real})
    predict_counterfactual(fit, Xpre)
end

"""
    residuals(fit, Xpre, Ypre)

Return pre-treatment residuals `Ypre - fitted_pre(fit, Xpre)`.
"""
function residuals(fit::MSCFit, Xpre::AbstractMatrix{<:Real}, Ypre::AbstractMatrix{<:Real})
    Yf = Matrix{Float64}(Ypre)
    fitted = fitted_pre(fit, Xpre)
    size(Yf) == size(fitted) || throw(DimensionMismatch("Ypre and fitted values must have the same size"))
    return Yf - fitted
end

"""
    prefit_rmse(fit, Xpre, Ypre)

Root mean squared pre-treatment fitting error.
"""
prefit_rmse(fit::MSCFit, Xpre::AbstractMatrix{<:Real}, Ypre::AbstractMatrix{<:Real}) =
    sqrt(mean(abs2, residuals(fit, Xpre, Ypre)))

function _msrl_single(X::Matrix{Float64}, Y::Matrix{Float64}, lambda::Float64;
                      theta0::Union{Nothing,Matrix{Float64}}=nothing,
                      max_iter::Int=10_000, tol::Float64=1e-6,
                      backtrack::Float64=0.5, verbose::Bool=false)
    nobs, ncontrols = size(X)
    nobs == size(Y, 1) || throw(DimensionMismatch("X and Y must have the same number of rows"))
    ntreated = size(Y, 2)
    theta_prev = isnothing(theta0) ? zeros(ncontrols, ntreated) : copy(theta0)
    theta = copy(theta_prev)
    alpha_prev = 1.0
    alpha = 1.0
    lip = max(opnorm(X)^2 / sqrt(nobs), eps(Float64))
    obj_values = Float64[]
    converged = false

    for iter in 1:max_iter
        z = theta + ((alpha_prev - 1.0) / alpha) * (theta - theta_prev)
        resid = Y - X * z
        F = svd(resid; full=false)
        grad = -(X' * (F.U * F.Vt)) / sqrt(nobs)

        obj_z = nuclear_norm(resid) / sqrt(nobs) + lambda * sum(abs, z)
        local candidate
        local obj_candidate
        local step_lip = lip

        while true
            candidate = soft_threshold(z - grad / step_lip, lambda / step_lip)
            diff = candidate - z
            candidate_resid = Y - X * candidate
            obj_candidate = nuclear_norm(candidate_resid) / sqrt(nobs) + lambda * sum(abs, candidate)

            quad = obj_z + sum(grad .* diff) + 0.5 * step_lip * sum(abs2, diff) +
                   lambda * (sum(abs, candidate) - sum(abs, z))
            if obj_candidate <= quad + 1e-10 || step_lip > 1e16
                break
            end
            step_lip /= backtrack
        end

        push!(obj_values, obj_candidate)
        if verbose && (iter == 1 || iter % 100 == 0)
            @info "MSC iteration" iter obj=obj_candidate nonzero=count(!iszero, candidate)
        end

        denom = max(1.0, norm(theta))
        if norm(candidate - theta) / denom <= tol
            theta_prev = theta
            theta = candidate
            converged = true
            return theta, obj_values, iter, converged
        end

        theta_prev = theta
        theta = candidate
        alpha_prev = alpha
        alpha = (1.0 + sqrt(1.0 + 4.0 * alpha^2)) / 2.0
        lip = step_lip
    end

    return theta, obj_values, max_iter, converged
end

function _fit_path(X::Matrix{Float64}, Y::Matrix{Float64}, lambdas::AbstractVector{<:Real};
                   max_iter::Int, tol::Float64, verbose::Bool)
    theta = zeros(size(X, 2), size(Y, 2))
    fits = Matrix{Float64}[]
    obj_path = Vector{Float64}[]
    iters = Int[]
    conv = Bool[]
    for λ in lambdas
        theta, obj, iter, converged = _msrl_single(
            X, Y, Float64(λ);
            theta0=theta,
            max_iter=max_iter,
            tol=tol,
            verbose=verbose,
        )
        push!(fits, copy(theta))
        push!(obj_path, obj)
        push!(iters, iter)
        push!(conv, converged)
    end
    return fits, obj_path, iters, conv
end

function _folds(nobs::Int, nfolds::Int, rng::AbstractRNG)
    2 <= nfolds <= nobs || throw(ArgumentError("nfolds must be between 2 and the number of rows"))
    perm = randperm(rng, nobs)
    folds = [Int[] for _ in 1:nfolds]
    for (i, row) in enumerate(perm)
        push!(folds[mod1(i, nfolds)], row)
    end
    return folds
end

function _cv_errors(X::Matrix{Float64}, Y::Matrix{Float64}, lambdas::AbstractVector{<:Real};
                    nfolds::Int, rng::AbstractRNG, max_iter::Int, tol::Float64)
    folds = _folds(size(X, 1), nfolds, rng)
    errors = zeros(length(lambdas))
    for fold in folds
        test = sort(fold)
        train_mask = trues(size(X, 1))
        train_mask[test] .= false
        Xtrain, Ytrain = X[train_mask, :], Y[train_mask, :]
        Xtest, Ytest = X[test, :], Y[test, :]
        thetas, _, _, _ = _fit_path(Xtrain, Ytrain, lambdas; max_iter=max_iter, tol=tol, verbose=false)
        for (j, θ) in enumerate(thetas)
            errors[j] += sum(abs2, Ytest - Xtest * θ)
        end
    end
    return errors ./ nfolds
end

function _as_dataframe(table)
    table isa DataFrame && return copy(table)
    Tables.istable(table) || throw(ArgumentError("input must be a DataFrame or Tables.jl-compatible table"))
    return DataFrame(table)
end

function _check_columns(df::DataFrame, cols::Vector{Symbol})
    missing_cols = setdiff(cols, Symbol.(names(df)))
    isempty(missing_cols) || throw(ArgumentError("missing columns: $(missing_cols)"))
end

"""
    panel_matrices(table, unit, time, outcome, treatment; treated_last=true)

Convert a balanced long panel to `MSCPanel`. Treatment must follow a common-adoption
design: controls are never treated, treated units are untreated through `T0` and treated
in every post-treatment period.
"""
function panel_matrices(table, unit::Symbol, time::Symbol, outcome::Symbol, treatment::Symbol;
                        treated_last::Bool=true)
    df = _as_dataframe(table)
    _check_columns(df, [unit, time, outcome, treatment])
    df = select(df, [unit, time, outcome, treatment])
    if any(ismissing, eachcol(df) |> Iterators.flatten)
        throw(ArgumentError("panel_matrices does not allow missing values"))
    end

    units0 = sort(collect(unique(df[!, unit])))
    times = sort(collect(unique(df[!, time])))
    n = length(units0)
    t = length(times)
    nrow(df) == n * t || throw(ArgumentError("panel must be balanced with one row per unit-time cell"))

    sort!(df, [unit, time])
    Y0 = Matrix{Float64}(reshape(Float64.(df[!, outcome]), t, n)')
    W0 = Matrix{Float64}(reshape(Float64.(df[!, treatment]), t, n)')
    all(x -> x == 0.0 || x == 1.0, W0) || throw(ArgumentError("treatment must be binary 0/1"))

    treated = vec(any(W0 .== 1.0; dims=2))
    any(treated) || throw(ArgumentError("no treated units found"))
    all(W0[.!treated, :] .== 0.0) || throw(ArgumentError("control units must never be treated"))

    first_treated = findfirst(vec(any(W0 .== 1.0; dims=1)))
    isnothing(first_treated) && throw(ArgumentError("no post-treatment period found"))
    T0 = first_treated - 1
    T0 >= 1 || throw(ArgumentError("at least one pre-treatment period is required"))
    T0 < t || throw(ArgumentError("at least one post-treatment period is required"))
    all(W0[:, 1:T0] .== 0.0) || throw(ArgumentError("no unit may be treated before T0"))
    all(W0[treated, (T0 + 1):end] .== 1.0) ||
        throw(ArgumentError("treated units must stay treated after common adoption"))

    if treated_last
        order = vcat(findall(.!treated), findall(treated))
    else
        order = collect(1:n)
    end

    Y = Y0[order, :]
    W = W0[order, :]
    units = Any[units0[i] for i in order]
    N0 = sum(.!treated)
    return MSCPanel(Y, N0, T0, W, units, Any[t for t in times])
end

"""
    simulate_msc(; kwargs...)

Generate a low-dimensional synthetic-control simulation with many treated units.
Returns a named tuple with `Y`, `N0`, `T0`, true `theta`, untreated counterfactuals,
and metadata.
"""
function simulate_msc(; ncontrols::Int=30, ntreated::Int=8, T0::Int=60, T1::Int=10,
                      rank::Int=3, sparsity::Int=5, tau::Real=1.0,
                      noise::Real=0.05, seed::Union{Nothing,Int}=nothing)
    ncontrols > 0 || throw(ArgumentError("ncontrols must be positive"))
    ntreated > 0 || throw(ArgumentError("ntreated must be positive"))
    T0 > 1 || throw(ArgumentError("T0 must be greater than 1"))
    T1 > 0 || throw(ArgumentError("T1 must be positive"))
    1 <= sparsity <= ncontrols || throw(ArgumentError("sparsity must be between 1 and ncontrols"))
    rng = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)
    T = T0 + T1

    factors = randn(rng, T, rank)
    loadings = randn(rng, rank, ncontrols)
    X = factors * loadings .+ Float64(noise) .* randn(rng, T, ncontrols)

    theta = zeros(ncontrols, ntreated)
    for j in 1:ntreated
        donors = randperm(rng, ncontrols)[1:sparsity]
        raw = rand(rng, sparsity)
        signs = rand(rng, [-1.0, 1.0], sparsity)
        theta[donors, j] = signs .* raw ./ sum(raw)
    end

    untreated = X * theta .+ Float64(noise) .* randn(rng, T, ntreated)
    treated = copy(untreated)
    treated[(T0 + 1):end, :] .+= Float64(tau)
    Y = vcat(X', treated')
    W = vcat(zeros(ncontrols, T), hcat(zeros(ntreated, T0), ones(ntreated, T1)))
    units = Any[vcat(["control_$i" for i in 1:ncontrols], ["treated_$i" for i in 1:ntreated])...]
    times = Any[collect(1:T)...]
    return (
        Y=Y,
        N0=ncontrols,
        T0=T0,
        W=W,
        theta=theta,
        untreated=untreated,
        treated=treated,
        tau=Float64(tau),
        units=units,
        times=times,
    )
end

"""
    fit_msc(Xpre, Ypre; kwargs...)

Fit the multivariate synthetic-control model from pre-treatment matrices.

`Xpre` is `T0 x n_controls`; `Ypre` is `T0 x n_treated`. The returned
coefficient matrix maps controls to treated counterfactuals.
"""
function fit_msc(Xpre::AbstractMatrix{<:Real}, Ypre::AbstractMatrix{<:Real};
                 lambdas::Union{Nothing,AbstractVector{<:Real}}=nothing,
                 nlambda::Int=50,
                 lambda_min_ratio::Float64=1e-2,
                 nfolds::Union{Nothing,Int}=nothing,
                 standardize::Bool=false,
                 max_iter::Int=10_000,
                 tol::Float64=1e-6,
                 rng::AbstractRNG=Random.default_rng(),
                 verbose::Bool=false)
    size(Xpre, 1) == size(Ypre, 1) ||
        throw(DimensionMismatch("Xpre and Ypre must have the same number of rows"))
    size(Xpre, 1) >= 2 || throw(ArgumentError("at least two pre-treatment periods are required"))

    X, Y, x_mean, x_scale, y_mean = _center_scale(Xpre, Ypre, standardize)
    λs = isnothing(lambdas) ? _lambda_path(X, Y; nlambda=nlambda, lambda_min_ratio=lambda_min_ratio) :
         Float64.(collect(lambdas))
    all(λ -> λ >= 0 && isfinite(λ), λs) || throw(ArgumentError("lambdas must be finite and nonnegative"))
    isempty(λs) && throw(ArgumentError("lambdas cannot be empty"))

    cv = isnothing(nfolds) ? nothing :
         _cv_errors(X, Y, λs; nfolds=nfolds, rng=rng, max_iter=max_iter, tol=tol)
    best = isnothing(cv) ? length(λs) : argmin(cv)

    thetas, obj_paths, iters, conv = _fit_path(X, Y, λs; max_iter=max_iter, tol=tol, verbose=verbose)
    theta_scaled = thetas[best]
    theta = theta_scaled ./ reshape(x_scale, :, 1)
    obj = obj_paths[best]

    return MSCFit(
        theta,
        λs[best],
        λs,
        cv,
        obj,
        iters[best],
        conv[best],
        x_mean,
        x_scale,
        y_mean,
        standardize,
    )
end

"""
    predict_counterfactual(fit, Xpost)

Predict untreated counterfactual outcomes for treated units. `Xpost` is
`Tpost x n_controls`; the result is `Tpost x n_treated`.
"""
function predict_counterfactual(fit::MSCFit, Xpost::AbstractMatrix{<:Real})
    Xf = Matrix{Float64}(Xpost)
    size(Xf, 2) == size(fit.theta, 1) ||
        throw(DimensionMismatch("Xpost must have one column per control unit"))
    intercept = fit.y_mean - vec(fit.x_mean' * fit.theta)
    return Xf * fit.theta .+ reshape(intercept, 1, :)
end

"""
    att(fit, Xpost, Ytreated_post)

Estimate ATT from a fitted MSC model and post-treatment outcomes.
"""
function att(fit::MSCFit, Xpost::AbstractMatrix{<:Real}, Ytreated_post::AbstractMatrix{<:Real})
    y0 = predict_counterfactual(fit, Xpost)
    Yf = Matrix{Float64}(Ytreated_post)
    size(Yf) == size(y0) || throw(DimensionMismatch("Ytreated_post must match predicted counterfactual size"))
    return mean(Yf .- y0)
end

"""
    effect_matrix(est)

Return post-treatment effects by period and treated unit.
"""
effect_matrix(est::MSCEstimate) = est.treated_post - est.counterfactual

"""
    effect_curve(est)

Average treatment effect by post-treatment period.
"""
effect_curve(est::MSCEstimate) = vec(mean(effect_matrix(est); dims=2))

"""
    unit_effects(est)

Average treatment effect by treated unit.
"""
unit_effects(est::MSCEstimate) = vec(mean(effect_matrix(est); dims=1))

"""
    weights_table(fit; units=nothing, treated_units=nothing)

Return the estimated weight matrix as a DataFrame with one row per nonzero
control-treated pair.
"""
function weights_table(fit::MSCFit; units=nothing, treated_units=nothing, atol::Real=0.0)
    θ = fit.theta
    ncontrols, ntreated = size(θ)
    ctrl = isnothing(units) ? Any["control_$i" for i in 1:ncontrols] : Any[units...]
    trt = isnothing(treated_units) ? Any["treated_$j" for j in 1:ntreated] : Any[treated_units...]
    length(ctrl) == ncontrols || throw(ArgumentError("units must have one label per control"))
    length(trt) == ntreated || throw(ArgumentError("treated_units must have one label per treated unit"))

    control_col = Any[]
    treated_col = Any[]
    weight_col = Float64[]
    for j in 1:ntreated, i in 1:ncontrols
        if abs(θ[i, j]) > atol
            push!(control_col, ctrl[i])
            push!(treated_col, trt[j])
            push!(weight_col, θ[i, j])
        end
    end
    return DataFrame(control=control_col, treated=treated_col, weight=weight_col)
end

weights_table(est::MSCEstimate; atol::Real=0.0) = weights_table(
    est.fit;
    units=est.units[1:est.N0],
    treated_units=est.units[(est.N0 + 1):end],
    atol=atol,
)

"""
    cv_table(fit)

Return cross-validation errors by lambda. If the model was fit without CV, the
error column is `missing`.
"""
function cv_table(fit::MSCFit)
    err = isnothing(fit.cv_errors) ? fill(missing, length(fit.lambdas)) : fit.cv_errors
    selected = fit.lambdas .== fit.lambda
    return DataFrame(lambda=fit.lambdas, cv_error=err, selected=selected)
end

"""
    msc_estimate(Y, N0, T0; kwargs...)

Estimate ATT from an outcome panel matrix. Rows `1:N0` are controls, remaining rows
are treated. Columns `1:T0` are pre-treatment periods.
"""
function msc_estimate(Y::AbstractMatrix{<:Real}, N0::Integer, T0::Integer; kwargs...)
    Yf = Matrix{Float64}(Y)
    n, t = size(Yf)
    1 <= N0 < n || throw(ArgumentError("N0 must leave at least one control and one treated unit"))
    1 <= T0 < t || throw(ArgumentError("T0 must leave at least one pre- and post-treatment period"))

    Xpre = Yf[1:N0, 1:T0]'
    Ypre = Yf[(N0 + 1):end, 1:T0]'
    fit = fit_msc(Xpre, Ypre; kwargs...)

    Xpost = Yf[1:N0, (T0 + 1):end]'
    Ytreated_post = Yf[(N0 + 1):end, (T0 + 1):end]'
    y0 = predict_counterfactual(fit, Xpost)
    units = Any[vcat(["control_$i" for i in 1:N0], ["treated_$i" for i in 1:(n - N0)])...]
    times = Any[collect(1:t)...]
    return MSCEstimate(mean(Ytreated_post .- y0), fit, y0, Ytreated_post, Yf, Int(N0), Int(T0), units, times)
end

function msc_estimate(panel::MSCPanel; kwargs...)
    est = msc_estimate(panel.Y, panel.N0, panel.T0; kwargs...)
    return MSCEstimate(
        est.estimate,
        est.fit,
        est.counterfactual,
        est.treated_post,
        panel.Y,
        panel.N0,
        panel.T0,
        panel.units,
        panel.times,
    )
end

function msc_estimate(table, unit::Symbol, time::Symbol, outcome::Symbol, treatment::Symbol; kwargs...)
    panel = panel_matrices(table, unit, time, outcome, treatment)
    return msc_estimate(panel; kwargs...)
end

function _subset_panel(Y::Matrix{Float64}, N0::Int, control_idx::Vector{Int}, treated_idx::Vector{Int})
    controls = control_idx
    treated = treated_idx
    return vcat(Y[controls, :], Y[treated, :]), length(controls)
end

"""
    placebo(est; T0_placebo=nothing, kwargs...)

Move the treatment date into the pre-treatment period and estimate a placebo ATT.
"""
function placebo(est::MSCEstimate; T0_placebo::Union{Nothing,Int}=nothing, kwargs...)
    placebo(est.Y, est.N0, est.T0; T0_placebo=T0_placebo, kwargs...)
end

"""
    placebo(Y, N0, T0; T0_placebo=floor(Int, T0/2), kwargs...)

Estimate a time-placebo effect using only the original pre-treatment periods.
"""
function placebo(Y::AbstractMatrix{<:Real}, N0::Integer, T0::Integer;
                 T0_placebo::Union{Nothing,Int}=nothing, kwargs...)
    Yf = Matrix{Float64}(Y)
    T0p = isnothing(T0_placebo) ? max(1, floor(Int, T0 / 2)) : T0_placebo
    1 <= T0p < T0 || throw(ArgumentError("T0_placebo must be between 1 and T0 - 1"))
    return msc_estimate(Yf[:, 1:T0], N0, T0p; kwargs...)
end

"""
    placebo_distribution(Y, N0, T0; replications=100, ntreated=nothing, kwargs...)

Randomly assign control units to placebo treatment groups and estimate MSC effects
using the original treatment date. This is most useful when `N0` is much larger than
the number of treated units.
"""
function placebo_distribution(Y::AbstractMatrix{<:Real}, N0::Integer, T0::Integer;
                              replications::Int=100,
                              ntreated::Union{Nothing,Int}=nothing,
                              rng::AbstractRNG=Random.default_rng(),
                              kwargs...)
    Yf = Matrix{Float64}(Y)
    n, _ = size(Yf)
    N1 = n - N0
    ntrt = isnothing(ntreated) ? N1 : ntreated
    N0 > ntrt || throw(ArgumentError("need more controls than placebo treated units"))
    replications > 0 || throw(ArgumentError("replications must be positive"))
    estimates = zeros(replications)
    for r in 1:replications
        perm = randperm(rng, N0)
        placebo_treated = sort(perm[1:ntrt])
        placebo_controls = sort(perm[(ntrt + 1):end])
        Yp, N0p = _subset_panel(Yf, N0, placebo_controls, placebo_treated)
        estimates[r] = msc_estimate(Yp, N0p, T0; kwargs...).estimate
    end
    return estimates
end

placebo_distribution(est::MSCEstimate; kwargs...) =
    placebo_distribution(est.Y, est.N0, est.T0; kwargs...)

"""
    placebo_se(Y, N0, T0; kwargs...)

Standard deviation of `placebo_distribution`.
"""
placebo_se(Y::AbstractMatrix{<:Real}, N0::Integer, T0::Integer; kwargs...) =
    std(placebo_distribution(Y, N0, T0; kwargs...))

placebo_se(est::MSCEstimate; kwargs...) = placebo_se(est.Y, est.N0, est.T0; kwargs...)

@recipe function f(est::MSCEstimate)
    tpost = est.times[(length(est.times) - size(est.treated_post, 1) + 1):end]
    eff = effect_curve(est)
    xlabel --> "Time"
    ylabel --> "Effect"
    label --> "Average effect"
    seriestype --> :line
    linewidth --> 2
    tpost, eff
end

@recipe function f(fit::MSCFit)
    if isnothing(fit.cv_errors)
        x = collect(1:length(fit.objective_values))
        y = fit.objective_values
        xlabel --> "Iteration"
        ylabel --> "Objective"
        label --> "Objective"
        seriestype --> :line
        return x, y
    end
    tbl = cv_table(fit)
    xlabel --> "lambda"
    ylabel --> "CV error"
    xscale --> :log10
    label --> "CV error"
    seriestype --> :line
    tbl.lambda, Float64.(tbl.cv_error)
end

end
