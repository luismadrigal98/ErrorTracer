# ErrorTracer 1.2.0 (development)

This release addresses the correctness issues raised in peer review, starting
with the uncertainty-budget consistency of the decomposition.

## Breaking changes / behaviour changes

* **Total variance now contains environmental uncertainty (T1).**
  `decompose_uncertainty()` previously computed `total_var` from
  `posterior_predict()` on the *unperturbed* predictors, so `env_var` was an
  add-on that sat *outside* the total — the reported components could sum to
  well over 100% of `total_var`. `total_var` is now **defined as the sum of
  its components** (law of total variance):
  `param_var + env_var + residual_var` for iid models, plus `temporal_var`
  for autocorrelation models. The percentage shares now sum to 100% exactly,
  and `env_var` is a genuine sub-share of a total that contains it.

* **`et_predict()` folds environmental uncertainty into the credible interval
  by default when `env_noise` is supplied.** The `include_env_in_ci` default
  changed from `FALSE` to `NULL` (auto): env is now included whenever the
  caller passes non-zero `env_noise`, so the reported forecast interval —
  and therefore `shelf_life()` and `et_calibrate()` coverage — reflect
  predictor/driver uncertainty and stay coherent with `total_var`. Pass
  `include_env_in_ci = FALSE` to recover the previous parameter+residual
  interval. The env-inclusive predictive is now constructed family-generally
  (re-centring each posterior-predictive residual on the perturbed mean),
  rather than the previous Gaussian-only `lp + sigma * N(0,1)`.

* `et_sensitivity_profile()` / `et_plot_sensitivity()`: `env_share` is now
  `env_var / total_var` (previously `env_var / (env_var + total_var)`), which
  double-counted env once env entered the total.

* **`extract_priors()` gains a `shrinkage` argument, defaulting to `"zero"`
  (T5).** This addresses the prior/likelihood double-use objection at its
  source. The previous behaviour — using the fitted coefficients as prior
  *means* and then refitting the Bayesian model on the same data — reused the
  data twice and could make posteriors overconfident. The new default centres
  every coefficient prior at 0 (a regularizing, ridge-style prior) with a width
  scaled by the coefficient *magnitude* `multiplier * |coef|` (floored at
  `min_sd`), so a well-estimated large coefficient is not shrunk toward 0 while
  its estimated *location* is no longer reused. Set `shrinkage = "estimate"` to
  recover the previous informative-mean behaviour (with a documented double-use
  caveat). For `ranger` the prior mean is always 0, so `shrinkage` has no
  effect there. The `et_prior_spec` object and its `print()` method now report
  the shrinkage mode.

## Bug fixes

* **AR/MA forecast variance now accumulates with lead time (T3).** For
  autocorrelation models (`ar()`/`ma()`/`arma()`/`cosy()`/`unstr()`/`sar()`/
  `car()`), `et_predict()` previously called `posterior_predict()` on the
  forecast rows alone, which does **not** propagate the residual
  autocorrelation across the horizon — the predictive variance (and hence the
  credible interval, `temporal_var`, and shelf life) stayed flat with lead
  time. `et_predict()` now continues the training series with the forecast
  rows (response blanked to `NA`) so brms iterates the fitted process forward:
  the forecast variance grows with lead time, integrates over the posterior of
  the AR coefficient (draws with |φ| ≥ 1 contribute super-linear growth,
  correctly handling the non-stationary / random-walk regime), and carries the
  initial-condition uncertainty from the forecast origin. The forecast rows'
  responses are always treated as unknown, so passing a hindcast's observed
  future values no longer leaks them into the interval.

* `et_predict()` caps `n_draws` at the number of draws in the fit (with a
  warning) instead of erroring, and now uses a shared, explicit set of draw
  ids across `posterior_predict()` and `posterior_linpred()` so the two are
  aligned draw-for-draw.

# ErrorTracer 1.1.0

## Breaking changes / deprecations

* `shelf_life()`: the `plausible_range` argument is deprecated. Use
  `response_scale` instead. `plausible_range` still works but emits a
  warning. It will be removed in a future release.

## New features

* **PIT / rank-histogram calibration diagnostics (T4).** New `et_pit()` computes
  the probability integral transform of held-out observations under the reported
  posterior predictive distribution (the same distribution the intervals and
  `et_calibrate()` coverage use, environmental uncertainty included), and
  `et_plot_pit()` renders the PIT histogram with a uniform reference band. This
  is the full-distribution diagnostic — a U-shape flags overconfident
  intervals, a hump flags under-confidence, a tilt flags bias — that coverage
  alone cannot provide. `et_pit()` supports randomized PIT for discrete
  predictions. `et_prediction` objects now store `predictive_draws`, the exact
  draws the interval / PIT / calibration all share.

* **Pathwise driver ensembles for genuine reforecasts (T2).** `et_predict()`
  gains an `env_ensemble` argument: a named list of \eqn{M \times H} matrices
  (M scenario trajectories over the H forecast steps), one per uncertain driver
  predictor. Each posterior draw is paired with a whole covariate *trajectory*,
  so driver uncertainty is temporally coherent and **accumulates with lead
  time** as a real driver ensemble (CMIP members, resampled climatology) fans
  out — unlike the existing `env_noise`, whose independent per-observation
  jitter averages out and leaves the driver variance flat. This is the
  recommended way to build a reforecast that uses only information available at
  the issue date, and it reproduces the growing driver-variance fraction of
  Dietze (2017) and Thomas et al. (2020). `env_ensemble` supersedes `env_noise`
  for the perturbation and is folded into the credible interval by default.

* `shelf_life()` now includes `se_t_star` in the projected-mode horizon
  attribute. This is the delta-method standard error of the projected
  crossing time `t* = (τ − a) / b`, propagating the linear-fit covariance.

* `decompose_uncertainty()` and the underlying `.decompose_from_arrays()`
  helper now support non-Gaussian families (Binomial, Poisson, Student-t,
  Negative-Binomial, Beta, Gamma) by computing all variance components on
  the response scale via the inverse link. For Gaussian identity the result
  is numerically unchanged.

* `et_predict()` gains an `n_env_draws` argument. Setting it > 1 averages
  multiple independent perturbations per posterior draw, reducing Monte
  Carlo noise in the environmental variance estimate. The decomposition
  data frame gains a `v_env_mcse` column reporting the chi-squared standard
  error of `env_var`.

* `decompose_uncertainty()` now reports a fourth `temporal_var` component
  when the model formula contains an autocorrelation term (`ar()`,
  `ma()`, `arma()`, `cosy()`, `unstr()`, `sar()`, or `car()`). The
  component is computed as `pmax(0, total_var - (param_var + env_var +
  residual_var))` and captures the autocorrelation-induced predictive
  variance beyond the iid sum, so the four components reconstruct
  `total_var` modulo Monte Carlo error. `residual_var` for autocor
  models is interpreted as the innovation variance (not the stationary
  marginal variance). `et_plot_decomposition()` adds a Temporal segment
  to the stacked bars automatically when the column is present, and
  `print.et_prediction()` includes the new row in its summary.

## Bug fixes

* `et_plot_calibration()` previously only recognised a column literally
  named `group` as the sub-group identifier. Calibration data frames built
  by binding per-group results with descriptive column names
  (e.g.\ `cluster_id`, `species`, `regime`) were silently collapsed into a
  single un-grouped series, producing plots with multiple overlapping
  points per nominal level and a zig-zagging connecting line. The function
  now auto-detects any single non-canonical column (anything other than
  `ci_level`, `nominal`, `observed_coverage`, `n_obs`, `calibration_error`,
  `sharpness`) with more than one unique value and uses it as the grouping
  variable. A new `group_col` argument allows the grouping column to be
  set explicitly, or `group_col = NA` to force a single un-grouped series.

# ErrorTracer 1.0.0

* Initial CRAN submission.

* Initial CRAN release.
* Added full Bayesian error propagation pipeline using `brms`.
* Added three-way variance decomposition and forecast shelf life metrics.
