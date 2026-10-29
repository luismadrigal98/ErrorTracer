## Summary

This release (1.4.2) adds one feature and fixes a portability bug. The version
on CRAN is 1.3.1, so this submission carries the 1.4.1 and 1.4.0 material as
well; neither was submitted. The 1.4.2 changes come first, then 1.4.1, then
1.4.0.

**Please note two changes in default output and one in reproducibility.**
`decompose_uncertainty()` now reports `total_var` as the *measured* posterior
predictive variance rather than the sum of its channels, and `temporal_var` as
a signed remainder rather than a floored one. `et_predict()` now seeds its own
draws. All three correct results that were wrong rather than merely
conservative: **any variance budget, credible interval, PIT, calibration or
shelf life produced by 1.3.x for a model carrying an autocorrelation term
should be recomputed rather than trusted.** Models without such a term are
affected only at the level of Monte Carlo error.

## 1.4.2 — a model fitted outside the package can now be used

`as_et_model()` wraps an existing `brms` fit in the object `et_fit()` returns,
and `et_predict()` / `et_diagnose()` accept a `brmsfit` directly, so the
prior-extraction and fitting stages are optional rather than a precondition for
the variance accounting. Bare Stan objects are refused with the
`brm(..., empty = TRUE)` route to take instead, since posterior draws alone
carry neither the family, the formula, nor the training data the decomposition
reads.

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
  than folding it into a channel, with `v_env_mcse` alongside `env_var` so the
  zero-floor on that channel warns when it fires beyond twice its Monte Carlo
  error.

## Test environments

* local Linux install (Ubuntu), R 4.6.1
* win-builder (devel and release)

## R CMD check results

`R CMD check --as-cran --run-donttest` on the built tarball:

0 errors | 0 warnings | 2 notes

**Note 1 — days since last update.** This flags submission cadence, not a code
issue. If it is too soon after the previous version, I am happy to wait and
resubmit; please tell me the preferred interval.

**Note 2 — slow examples.** The expected note for a Stan-backed package: five
examples each fit one Bayesian model, so almost all of that time is Stan model
**compilation** rather than computation. The examples are already minimal — 20
observations, one predictor, a single chain of 500 iterations, of which
sampling accounts for under a second — so the cost is not reducible by
shrinking them further. They are wrapped in `\donttest{}`. This note was
present in the previously accepted releases.

For context on total check time: the `testthat` suite is gated with
`skip_on_cran()`, so it runs in a few seconds under CRAN's settings (709
tests across 22 files locally with `NOT_CRAN=true`, 0 failures, 0 skips; the
suite grew from 572 at 1.3.1 as each defect below was pinned by a regression
test). Vignettes do not fit models at build time.

## Reverse dependencies

None on CRAN.
