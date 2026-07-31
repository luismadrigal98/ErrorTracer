# ErrorTracer 1.3.1

Corrects `shelf_life()`, which could report a forecast horizon that rested on a
single noisy period. **Any reported shelf life should be re-checked**; series
whose CI-width / response-scale ratio degrades monotonically are unaffected.

## Bug fixes

* **`shelf_life()` reported the first threshold exceedance as the horizon.**
  When the ratio hovers near `threshold` without trending, one period can
  exceed it by chance, and that isolated excursion was promoted to a "shelf
  life". Observed on a real series whose ratio sat at ~0.92 across a 16-period
  window with a single excursion to 1.038: the reported horizon was not
  reproducible across sampler settings, varying between "no crossing", a
  projection far outside the window, and two different crossing periods, from
  MCMC noise alone.

  `shelf_life()` now requires the exceedance to **persist**. The new `min_run`
  argument (default `2`) sets how many consecutive uninformative periods
  constitute a crossing; `min_run = 1` restores the previous behaviour. A
  crossing already in force at the first lead is still detected, and a dip
  inside an otherwise degrading series no longer resets the horizon.

* **The projection mode extrapolated slopes indistinguishable from flat.**
  `min_slope_for_projection` is a magnitude gate only, so a slope of ~1e-3 with
  a one-sided p-value of 0.42 cleared it and was extrapolated to a
  confident-looking crossing time. On the same real series above this produced a
  projected horizon 45 years out whose delta-method 95% interval spanned 328
  years and **included the past**. The projection now additionally requires the
  slope to be significantly positive, controlled by the new `projection_alpha`
  argument (default `0.05`); when it is not, the result is reported as a lower
  bound — the honest outcome for a forecast whose precision does not degrade
  in-window. Set `projection_alpha = 1` to restore the previous behaviour.

  The two fixes compose: a series that neither sustains a crossing nor trends
  detectably now returns `lower_bound` instead of either a spurious observed
  horizon or a spurious projection.

* **`et_skill_score()` had the same first-crossing defect in its forecast
  limit.** The null-relative forecast limit was the first lead time with
  `skill <= 0`, with no persistence requirement. Skill hovering near zero
  crosses by chance, so one unlucky lead was promoted to "the horizon" — found
  on a real comparison whose mean skill over the window was comfortably
  positive. `et_skill_score()` now takes `min_run` (default `2`, matching
  `shelf_life()`) and the returned `"forecast_limit"` attribute carries
  `n_below`, `first_below` and `min_run`, so isolated dips are visible rather
  than load-bearing. `min_run = 1` restores the previous behaviour.

## New features

* **Group-level variance channel for hierarchical models.** The decomposition
  now handles `(1 | g)`-style fits correctly, and distinguishes the two budgets
  that a hierarchical model implies:

  - Predicting for a group the fit has **seen**: the group effect is an
    estimated parameter, its posterior uncertainty is part of `param_var`, and
    there is no separate channel. Conditioning on a known group *reduces*
    predictive variance relative to the population level, so the naive
    `Var(full) - Var(population)` is negative and would be meaningless as a
    "group variance" (measured on a 6-group fixture: 0.022 vs 0.244).
  - Predicting for a **new** level: the budget integrates over the group-level
    distribution and a `group_var` column appears, carrying the between-group
    variance. Verified to recover the fitted `tau^2` on a known-truth fixture.

  `et_predict()` detects new levels automatically and reports which budget it is
  computing. `param_var + group_var + env_var + residual_var (+ temporal_var)`
  still equals `total_var` exactly.

* **`et_shelf_life_pool()` — pool shelf lives across forecast origins.** A
  horizon from a single forecast origin is a property of that training window as
  much as of the system. This pools horizons from repeated origins and returns a
  median with an interval. Origins where the forecast was still informative when
  the data ran out are **right-censored**, not missing: discarding them biases
  the pooled horizon downward, because the censored origins are exactly the ones
  with long horizons. Estimation is Kaplan–Meier via `survival` (a Recommended
  package, in `Suggests`), with a documented, deliberately interval-free
  fallback when it is unavailable. All origins censored is reported as a result
  — the forecast never became uninformative — not as an estimation failure.

* **`et_priors_split()` — remove the prior/likelihood double use entirely.**
  `shrinkage = "zero"` already stopped reusing point estimates as prior means,
  but the prior *scale*, and any variable *selection* performed by the
  first-stage model, still came from the rows that form the likelihood. This
  fits the first stage on a random (optionally stratified) subset and returns
  the disjoint complement for the Bayesian fit, so selection and likelihood
  share no observations and no post-selection caveat applies. The cost —
  a smaller likelihood sample — is documented rather than hidden. The
  user-supplied seed is restored afterwards so it cannot leak into the caller's
  RNG stream.

* **Exceedance diagnostics on every horizon.** The `horizon` attribute now
  carries `n_exceedances`, `frac_exceedance`, `first_exceedance` and `min_run`
  alongside `value` and `type`, so an isolated excursion is visible rather than
  silent. Comparing `first_exceedance` with `value` shows directly whether a
  horizon rests on a persistent crossing. Trend-mode horizons additionally carry
  `slope` and `slope_p`.

# ErrorTracer 1.3.0

Corrects the variance decomposition for models carrying a residual
autocorrelation term. **Any analysis using `ar()` / `ma()` / `arma()` /
`cosy()` / `unstr()` / `sar()` / `car()` must be re-run**; iid models are
unaffected (all three changes below are no-ops without an autocorrelation
term). Note that 1.2.1 ships the defect.


## Bug fixes — variance decomposition under autocorrelation

* **`env_var` was silently destroyed for any model with an `ar()` / `ma()` /
  `arma()` term.** Under brms's default `ar(cov = FALSE)` parameterisation the
  autocorrelation term enters the **mean** (`mu[n] += ar * err[n-1]`), not only
  the covariance. `decompose_uncertainty()` computed
  `env_var = Var(mu_perturbed) - Var(mu_draws)` where `mu_perturbed` came from an
  internal hand-rolled linear predictor (no AR term) while `mu_draws` came from
  `posterior_linpred()` (AR term included). The difference was therefore
  `V_env - V_AR`, which turns negative as the AR error grows and was then clamped
  to zero by `pmax(., 0)`.

  In a controlled test (AR(1), phi = 0.9, analytic `V_env = beta^2 * delta^2 =
  0.152`) the reported `env_var` averaged **0.012 and was floored to exactly zero
  at 13 of 15 lead times**; the iid control recovered 0.173 correctly. Both
  `posterior_linpred()` calls now pass `incl_autocor = FALSE`, which returns the
  pure regression linear predictor — verified bit-identical to the internal
  predictor, so the AR term cancels exactly in the subtraction. `env_var` now
  averages 0.145 with no floored rows.

* **`param_var` no longer absorbs the autocorrelation error, and `temporal_var`
  no longer collapses to zero.** The same contamination meant `param_var` grew
  with lead time even on non-trending new data (0.154 -> 0.976 over 15 steps, an
  8.9x overstatement of parameter uncertainty at the far lead), which in turn
  made `temporal_var = pmax(0, pp_var - (param_var + residual_var))` evaluate to
  **zero at every lead for a series with phi = 0.82**. After the fix `param_var`
  is flat (0.154 -> 0.147) and `temporal_var` grows monotonically 0 -> 0.637,
  correctly reporting zero at lead 1 where a single innovation adds nothing
  beyond `residual_var`.

  **This changes reported numbers for every autocorrelation model.** Analyses
  using `ar()`/`ma()`/`arma()`/`cosy()`/`unstr()`/`sar()`/`car()` must be re-run.

* **The `env_var` zero-floor is now instrumented.** `decompose_uncertainty()`
  warns when the floor fires on rows whose shortfall exceeds twice its own Monte
  Carlo SE — too systematic to be jitter. The silent floor is what hid the bug
  above through two review rounds. The test is statistical rather than a raw
  count, so small prediction sets (where a floored row genuinely is jitter) do
  not trip it.

## Autocorrelation priors — treating the cause rather than the tail

* **`et_fit()` gains `ar_prior`, defaulting to `"weakly_informative"`.** `brms`
  places a *flat, unbounded* prior on `class = "ar"` / `"ma"`
  (`get_prior()` reports `(flat)` with no `lb`/`ub`), so nothing keeps the
  residual process inside the stationary region. On the short series typical of
  ecological panels the posterior routinely straddles the unit root, and the
  k-step forecast variance `sigma^2 (1 - phi^(2k)) / (1 - phi^2)` then diverges.
  The new default `normal(0, 0.5)` places ~95% of the prior mass inside
  `(-1, 1)` without hard-bounding, so a genuinely non-stationary series can
  still say so. `ar_prior = "stationary"` truncates to `(-1, 1)`;
  `ar_prior = "flat"` restores the old behaviour. A prior the caller supplies
  for these classes is never overridden.

  On a near-unit-root fixture (n = 23, true random-walk residual) this reduces
  the ratio of raw to 1%-winsorized predictive variance from **1072x to 1.7x**
  — i.e. the heavy tail that motivated the winsorization below was largely an
  artefact of the unbounded prior, not a property of the data.

* **`et_fit()` warns when more than 5% of an autocorrelation posterior sits at
  `|value| >= 1`**, reporting the exact fraction and posterior mean. For AR(1)
  this is precisely the posterior probability of a non-stationary process. The
  message directs the user to read long-lead intervals, `temporal_var` and
  `shelf_life()` as "no usable horizon" rather than as calibrated numbers.

## Variance-estimator consistency

* **`et_predict()` gains `var_trim` (default `0`), replacing the always-on,
  undocumented 1% winsorization.** Two problems are fixed:
  1. The winsorization was applied to the posterior-predictive variance *only*,
     while `param_var` and the environmental pair stayed raw. `temporal_var`
     was therefore a difference between two different estimators, biased low by
     whatever the winsorization removed (~3.6% for a Gaussian at `trim = 0.01`)
     — which blinded it to genuine autocorrelation shares below that threshold.
     `var_trim` is now applied uniformly to **every** channel.
  2. Reporting a winsorized variance next to a raw quantile credible interval
     meant the budget and the interval described *different distributions*. With
     the `var_trim = 0` default they describe the same one.

  The tail diagnostic is still computed (against a fixed 1% reference,
  regardless of `var_trim`, since it is a statement about the model rather than
  the estimator) and now recommends `ar_prior = "stationary"` instead of
  implying that a robust summary makes a diverging forecast trustworthy.

  **This changes reported variance components** wherever the old winsorization
  was biting — i.e. any near-unit-root autocorrelation fit.

* Regression tests added in `tests/testthat/test-decompose-autocor.R`, covering
  `env_var` recovery under `ar()`, floor-rate, no decay with lead time, positive
  `temporal_var`, flat `param_var`, and an iid control. Reproduction scripts in
  `tests/manual/`. The previous suite missed all of this because its mock builds
  `mu_perturbed` as `lp + noise`, i.e. it assumed the two arrays were already on
  the same footing.

# ErrorTracer 1.2.1

* **Fix: platform-dependent `v_env_mcse` on aarch64 (macOS arm64).** This is
  the failure that caused the CRAN archival on 2026-06-08.
  `decompose_uncertainty()` detected the "no environmental noise" case by
  testing the variance difference `v_env_raw == 0` exactly, and it computed the
  two column variances with different `stats::var()` `na.rm` settings. Those
  settings select different C code paths in `cov.c`, which round identically on
  x86_64 but not on aarch64. A last-bit difference therefore flipped the
  reported Monte Carlo SE of `env_var` from `0` to a large value (0.20 in the
  CRAN log) when `env_noise = 0`. The two variances are now computed the same
  way, and the zero test uses a relative tolerance instead of exact equality,
  so the result no longer depends on bit-level floating-point reproducibility.
  `env_var` itself was already correct to within 1e-10; only its Monte Carlo SE
  was affected. Regression tests added.

# ErrorTracer 1.2.0

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

## Documentation & robustness

* **Scope honesty and a convergence nudge (T8).** `et_fit()` now documents its
  boundary plainly (a thin `brms` wrapper; the decomposition of hierarchical
  fits is provisional pending a group-variance term, and `grouping=` fits
  independent per-group models rather than one multilevel model), and emits a
  warning on high Rhat / divergent transitions pointing to `et_diagnose()`
  before prediction. `extract_priors.ranger()` is documented and flagged
  (advisory) as \strong{experimental}: random-forest importance has no sign or
  GLM-coefficient scale, so its importance→prior-SD mapping is heuristic —
  prefer `glmnet`/`lm`/`glm` for a principled prior.

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

* **The AR/temporal decomposition is now tail-robust and consistent with the
  reported interval.** Two follow-on fixes to the autocorrelation forecast
  above. (i) `et_predict()` built the credible-interval / PIT / calibration
  draws (`predictive_draws`) from only the first `n_perturb` posterior draws
  while the variance decomposition used all `n_draws`, so the interval and the
  decomposition could describe different samples; the env-inflated predictive
  now spans all `n_draws` draws (the env perturbation is recycled across them).
  (ii) `temporal_var` (and hence `total_var`) is computed from a **winsorized**
  posterior-predictive variance, so that when the AR coefficient's posterior has
  mass at \eqn{|\phi| \ge 1} a handful of explosive, non-stationary draws can no
  longer dominate the estimate. Previously a single near-unit-root draw could
  inflate the raw variance by an order of magnitude and make `temporal_var`
  disagree wildly with the (quantile-based, tail-robust) credible interval
  reported alongside it. The decomposition now tracks that interval; a genuinely
  non-stationary AR posterior still yields a large temporal share. `et_predict()`
  emits a heads-up when the raw predictive variance greatly exceeds its robust
  value (the near-unit-root regime), pointing to `et_diagnose()`.

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

* **Optional global (Sobol) variance decomposition (T7).** New `et_sobol()`
  provides an order-independent alternative to the additive budget: a
  variance-based decomposition of the predictive mean into first-order
  \strong{parameter} and \strong{environmental} Sobol indices plus their
  \strong{interaction} (Saltelli 2010 / Jansen estimators, with output
  centring for Monte Carlo stability). The fast additive budget remains the
  default; `et_sobol()` is for the general case where a non-linear link or
  correlated drivers make the parameter\eqn{\times}driver interaction
  non-negligible, which the sequential additive split folds into its terms.

* **Null-model forecast skill + shelf-life skill gate (T6).** New
  `et_skill_score()` scores a forecast against a null model — a random walk
  (persistence) or climatology — with the continuous ranked probability score
  (CRPS) and reports the per-lead-time skill and the null-relative
  \strong{forecast limit} (the lead time at which the model stops beating the
  null), following the ecological-forecasting convention (Petchey et al. 2015;
  Wesselkamp et al. 2025; the NEON challenge). `shelf_life()` gains a `skill`
  argument: passing an `et_skill_score()` table gates the shelf life so a
  period counts as informative only if it is both precise \emph{and} skillful,
  preventing a precise-but-biased forecast from passing. `shelf_life()` is now
  documented as a \emph{precision} complement to the accuracy-relative forecast
  limit, and its `response_scale` now defaults to the training-response range
  (a documented, reproducible default) instead of requiring an arbitrary
  number.

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
