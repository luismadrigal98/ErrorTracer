## Summary

This is a bug-fix and feature release (1.3.0) following 1.2.1.

It corrects a substantive statistical error in the package's core
functionality: the variance decomposition returned incorrect components for
any model carrying a residual autocorrelation term. Users of `ar()` / `ma()` /
`arma()` / `cosy()` / `unstr()` / `sar()` / `car()` models should re-run their
analyses. Models without an autocorrelation term are unaffected — all three
changes below are provably no-ops in that case.

## The bug, and how it was found

`decompose_uncertainty()` computed the environmental variance component as
`Var(perturbed linear predictor) - Var(unperturbed linear predictor)`. The
first array was built internally from the posterior draws matrix and carried
no autocorrelation term; the second came from `brms::posterior_linpred()`,
which under brms's default `ar(cov = FALSE)` parameterisation *does* include
the autocorrelation contribution to the mean. The subtraction therefore
evaluated `V_env - V_AR`, which turns negative as the AR error grows and was
then clamped to zero by a `pmax(., 0)` guard intended only to absorb Monte
Carlo jitter.

The defect was invisible on inspection — the source comment asserted that
"the mean structure carries no autocorrelation" — and was found by testing the
package against a closed-form target instead. In a controlled AR(1) fit
(phi = 0.9) whose analytic environmental variance is
`beta^2 * delta^2 = 0.152`, the reported `env_var` averaged **0.012 and was
floored to exactly zero at 13 of 15 lead times**; the corresponding iid
control recovered 0.173 correctly. The same contamination inflated
`param_var` by a factor of 8.9 at the far lead and drove `temporal_var` to
**zero at every lead time** for a series with strong, unambiguous residual
autocorrelation.

## Fixes

* Both `posterior_linpred()` calls now pass `incl_autocor = FALSE`, returning
  the pure regression linear predictor. This is verified bit-identical to the
  internally-constructed predictor, so the autocorrelation term cancels
  exactly in the subtraction. `env_var` now recovers 0.145 against the 0.152
  target with no floored rows; `param_var` is flat in lead time as it should
  be; `temporal_var` grows monotonically and is correctly zero at lead 1.

* `var_trim` (new argument to `et_predict()`, default `0`) applies a single
  variance estimator to every channel. Previously the posterior-predictive
  variance was winsorized while the parameter and environmental terms were
  not, so `temporal_var` was a difference between two different estimators,
  biased low by roughly the amount the winsorization removed.

* `ar_prior` (new argument to `et_fit()`, default `"weakly_informative"`)
  replaces brms's flat, unbounded prior on autocorrelation coefficients with
  `normal(0, 0.5)`, and `et_fit()` now warns when more than 5% of the
  posterior lies at `|phi| >= 1`. On a short series the unbounded default
  routinely yields posteriors straddling the unit root, where the k-step
  forecast variance diverges. `ar_prior = "stationary"` truncates to
  `(-1, 1)`; `ar_prior = "flat"` restores the previous behaviour.

* `env_var`'s zero-floor is now instrumented: `decompose_uncertainty()` warns
  when the floor fires on rows whose shortfall exceeds twice its own Monte
  Carlo standard error, i.e. more systematically than jitter explains. The
  silent floor is what concealed the bug above.

New regression tests in `tests/testthat/test-decompose-autocor.R` fail on the
old behaviour and pass on the new. `tests/manual/` holds the two standalone
reproduction scripts. See NEWS.md for the full entry.

## Test environments

* local Linux install (Ubuntu), R 4.6.1
* win-builder (devel and release)

## R CMD check results

`R CMD check --as-cran --run-donttest` on the built tarball:

0 errors | 0 warnings | 1 note

The note is the expected one for a Stan-backed package:

```
* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
                        user system elapsed
decompose_uncertainty 55.095  2.676  58.800
et_calibrate          51.304  2.301  53.721
shelf_life            51.136  2.310  53.537
et_fit                49.957  2.267  52.356
```

These four examples each fit one Bayesian model, so almost all of that time is
Stan model **compilation** rather than computation. The examples are already
minimal — 20 observations, one predictor, a single chain of 500 iterations, of
which sampling accounts for under a second — so the cost is not reducible by
shrinking them further. They are wrapped in `\donttest{}`. This note was
present in the previously accepted 1.2.x releases and no example has been
added or enlarged since.

For context on total check time: the `testthat` suite is gated with
`skip_on_cran()`, so it runs in about 4 seconds under CRAN's settings (341
seconds locally with `NOT_CRAN=true`, as reported above). Vignettes do not fit
models at build time.

## Reverse dependencies

None on CRAN.
