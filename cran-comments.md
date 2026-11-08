## Summary

This release (1.4.3) changes one default so that the variance budget is as
tight as the fit allows. The version on CRAN is 1.3.1, so this submission also
carries the 1.4.2, 1.4.1 and 1.4.0 material. 1.4.1 and 1.4.0 were not
accepted; 1.4.2 was submitted on 5 November 2026 and is, as far as I know,
still awaiting manual inspection. If it is preferable for that submission to be
processed first, or for this one to replace it, I will do whichever the team
prefers. The changes are listed newest first.

**Please note two changes in default output and one in reproducibility.**
`decompose_uncertainty()` now reports `total_var` as the *measured* posterior
predictive variance rather than the sum of its components, and `temporal_var` as
a signed remainder rather than a floored one. `et_predict()` now seeds its own
draws. All three correct results that were wrong rather than merely
conservative: **any variance budget, credible interval, PIT, calibration or
shelf life produced by 1.3.x for a model carrying an autocorrelation term
should be recomputed rather than trusted.** Models without such a term are
affected only at the level of Monte Carlo error.

## 1.4.3 — the decomposition now uses every draw in the fit

`et_predict(n_draws =)` defaulted to 2000 draws even when the fit held many
more. Since 1.4.0 `total_var` is the *measured* variance of the predictive
draws rather than the sum of its components, so for an iid model the reported
`closure_gap` is exactly the Monte Carlo error of that measurement.
Subsampling therefore widened the gap for no benefit, because the draws are
paid for during sampling and prediction is matrix work. On a 90-observation
Gaussian fit the per-observation gap has mean zero throughout, but its spread
falls from 4.0% of `total_var` at 2000 draws to 1.2% at 8000.

The default is now `NULL`, meaning every available draw; an integer still
bounds memory on wide `newdata` and is capped at the draws present. **This
changes reported numbers**: a decomposition, interval, PIT, calibration or
shelf life from the old default is a noisier estimate of the same quantity,
not a different one.

### 1.4.3 also refuses a combination it used to accept

`brms` carries an autocorrelation term as response residuals for `gaussian` and
`student`, but as a vector of *latent* residuals for every other family under
its default `cov = FALSE`, and it cannot extend latent residuals to new rows.
So `et_fit(y ~ x + ar(time = t), family = poisson())` sampled successfully and
`et_predict()` then failed inside `brms` with "Cannot predict new latent
residuals". **`et_fit()` and `et_predict()` now stop on that combination**,
before Stan is invoked, naming `cov = TRUE` as the fix; `brms` offers the
covariance formulation at order one only, which the message also states.
Gaussian and student models are untouched: a Portal-style `ar(cov = FALSE)`
fit predicts bit-identically before and after this change. The `cov = TRUE`
path is permitted and now warns that its
driver-inclusive draws are not guaranteed to respect the response support: on a
small-count Poisson fit 12.3% were negative and all were non-integer.

## 1.4.2 — a model fitted outside the package can now be used

`as_et_model()` wraps an existing `brms` fit in the object `et_fit()` returns,
and `et_predict()` / `et_diagnose()` accept a `brmsfit` directly, so the
prior-extraction and fitting stages are optional rather than a precondition for
the variance accounting. Bare Stan objects are refused with the
`brm(..., empty = TRUE)` route to take instead, since posterior draws alone
carry neither the family, the formula, nor the training data the decomposition
reads.

Three changes to plot output. `et_plot_decomposition(proportional = TRUE)` now normalises each component by the components actually drawn, so the bars are true shares summing to 1 and the axis is a percentage; previously it divided by `total_var` while drawing only the named components, which left iid bars short of 1 by the closure gap. Tick labels are now angled based on the labels themselves rather than on the x variable, so long group names no longer collide.

Two further changes to plot output. `et_plot_*()` no longer sets a default title: these figures are made for manuscripts, where the caption is typeset by the journal, and a title in the image duplicates it. Axis labels, legend titles and facet strips are unchanged, and a title is still available through `+ labs(title = )`.

It also fixes a faceting defect in `et_plot_pit()`: the dashed uniform line
and the 95% consistency band were computed once from the mean group size and
drawn identically in every panel, so panels of unequal size were judged against
a reference that belonged to neither. Each panel now gets its own. No returned
value changes; the affected output is the plot.

One portability fix: `%||%` was used in `R/fit.R` and `R/predict.R` but never
defined or imported, so those calls resolved only on R >= 4.4.0, where base R
supplies it. `DESCRIPTION` declares R (>= 4.1.0), so the package would have
failed on the oldrel CRAN checks. It is now defined in `R/utils.R`.

## 1.4.1 — a prediction was not reproducible, and a plot axis was mislabelled

* **`et_predict()` drew from the ambient random stream.** The
  posterior-predictive draw and the driver perturbation both consume R's RNG
  and nothing seeded them: `et_fit(seed =)` reached Stan and stopped there. Two
  runs of the same fit on the same new data therefore differed, and by how much
  depended on what had used the stream earlier in the session. In this
  package's own case studies, two scripts that fit the same model to the same
  training window disagreed about one series' forecast by a full lead, because
  one reached `et_predict()` after a single fit and the other after several.
  `et_predict()` now takes `seed`, defaulting to the seed the model was fitted
  with, so a prediction is reproducible from the model alone; `seed = NA` keeps
  the ambient stream. The caller's stream is saved and restored, so seeding
  here never reruns anyone else's draws.

* **`et_plot_decomposition(group_col = )` plotted a sum of shares on an axis
  labelled as a share.** Each component was divided by its own row's
  `total_var` and the rows were then stacked at one x position per group, so a
  group with `H` observations produced a bar of height `H`. The colour
  proportions were unaffected — the picture was right and the axis was not.
  Rows are now aggregated before normalising, which also makes the bar equal
  the variance-weighted share that the package's own summaries report.

* **`et_plot_shelf_life()`** gained `max_ratio`, applied through
  `coord_cartesian()` so it zooms rather than filters, and its threshold label
  is now placed by clearance rather than pinned to the panel's left edge.

## 1.4.0 — three defects in the predictive draws

* **The predictor perturbation was recycled across posterior draws.**
  `et_predict()` perturbed only the first `n_perturb` draws and the inflation
  step then recycled that block across all `n_draws` rows, pairing the
  perturbed mean of one draw with the predictive residual of another. That
  destroys `Cov(mu, y - mu)`, which is approximately zero for an iid fit but
  substantially negative for an `ar()` fit, so the predictive variance was
  biased upward. On an annual ecological series the gap between predictive and
  posterior-predictive variance, which should equal `env_var`, came out at 20x,
  22x and 196x `env_var` for the three most persistent cases. Because those
  draws are the single source of truth for the credible interval, PIT,
  calibration *and* shelf life, the inflation reached all four: 90% intervals
  were 1-12% too wide. `n_perturb` now defaults to `n_draws`.

* **`binomial()` with `trials(n > 1)` decomposed on the probability scale while
  the predictive draws were counts.** With `n = 10`, `residual_var` came back
  ten times too small and `param_var` a hundred times too small, and the
  env-inflated draws added a probability to a count residual, landing outside
  the support. The predictive variance was effectively flat in the perturbation
  SD, so driver uncertainty was barely propagated at all.

* **The env-inclusive predictive draw was always a location shift.** It is now
  drawn from the response distribution at the perturbed mean, which coincides
  with the old construction under a Gaussian identity link and otherwise
  returns values on the family's own support. Models carrying a residual
  autocorrelation term keep the location shift deliberately, since the
  deviation is what carries the accumulation forward; the returned object
  records which construction was used.

## 1.4.0 — changes to the decomposition

* **`total_var` is measured, not defined**: it is `Var[posterior_predict]`, the
  spread of the distribution the reported interval is drawn from. A budget
  whose total is defined as the sum of its parts cannot report when the parts
  fail to account for it.
* **`temporal_var` is signed.** It is a covariance term, and a negative value
  is meaningful, since autocorrelated errors can reduce predictive variance
  relative to an independent-innovation baseline.
* **New `closure_gap` column for iid models**, reporting the shortfall rather
  than folding it into a component, with `v_env_mcse` alongside `env_var` so the
  zero-floor on that component warns when it fires beyond twice its Monte Carlo
  error.

## Test environments

* local Linux install (Ubuntu), R 4.6.2
* win-builder (devel and release)

## R CMD check results

`R CMD check --as-cran --run-donttest` on the built tarball, R 4.6.2 on Linux:

0 errors | 0 warnings | 2 notes

**Note, slow examples.** The expected note for a Stan-backed package: six
examples each fit one Bayesian model, so almost all of that time is Stan model
**compilation** rather than computation. The examples are already minimal (20
observations, one predictor, a single chain of 500 iterations, of which
sampling accounts for under a second), so the cost is not reducible by
shrinking them further. They are wrapped in `\donttest{}`. Measured elapsed
times for this submission: `as_et_model` 53.1s, `decompose_uncertainty`
50.5s, `shelf_life` 46.5s, `et_calibrate` 46.5s, `et_priors_split` 45.6s,
`et_fit` 45.5s (290s in total). An equivalent note was present in the
previously accepted releases.

**Note, submission cadence.** The incoming check reports "Days since last
update: 3. Number of updates in past 6 months: 7." That reflects
submission cadence rather than a code issue; if it is too soon I am happy to
wait and resubmit, and would welcome being told the preferred interval.

For context on total check time: the `testthat` suite is gated with
`skip_on_cran()`, so it runs in a few seconds under CRAN's settings (758
tests across 23 files locally with `NOT_CRAN=true`, 0 failures, 0 skips; the
suite grew from 572 at 1.3.1 as each defect below was pinned by a regression
test). Vignettes do not fit models at build time.

## CRAN incoming NOTE: possibly misspelled words

The incoming check flags `Sobol` and `Kaplan` in the Description. Both are
surnames and both are spelled correctly:

* **Sobol** — Ilya M. Sobol', after whom Sobol indices (variance-based global
  sensitivity indices) are named; the package implements them in
  `et_sobol()`.
* **Kaplan** — Edward L. Kaplan, of the Kaplan-Meier estimator, which
  `et_shelf_life_pool()` uses to pool right-censored shelf lives across
  rolling forecast origins.

Neither is a package or software name, so neither is quoted.

## Reverse dependencies

None on CRAN.
