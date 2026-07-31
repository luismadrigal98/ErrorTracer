## Summary

This is a bug-fix release (1.3.1). It supersedes 1.3.0, which was prepared but
never submitted, so this submission carries the fixes from both. 1.3.0's fix is
summarised second below; the new 1.3.1 material comes first.

**Please note a deliberate change in default output.** `shelf_life()` now
requires a threshold crossing to persist before reporting it as a forecast
horizon, and no longer extrapolates a trend that is statistically
indistinguishable from flat. Both changes correct results that were wrong
rather than merely conservative, so a horizon computed with 1.2.x or 1.3.0
should be recomputed rather than trusted. The previous behaviour remains
reachable via `min_run = 1` and `projection_alpha = 1`. I have kept this as a
patch release because it is a defect fix and the old behaviour is fully
recoverable, but I am flagging the behaviour change explicitly rather than
leaving a reviewer to discover it.

## 1.3.1 — two defects in `shelf_life()`

Both have the same shape: a noisy quantity promoted to a confident-looking
number.

* **The first threshold exceedance was reported as the horizon.** When the
  CI-width / response-scale ratio hovers near the threshold without trending, a
  single period exceeds it by chance. On a real 16-period series whose ratio sat
  at ~0.92 throughout with exactly one excursion to 1.038, that excursion was
  reported as a 12-period horizon. It was not reproducible: holding data, seed
  and priors fixed and varying only sampler settings moved the answer between
  "no crossing", a projection far outside the window, and two different crossing
  periods. `shelf_life()` now requires `min_run` consecutive uninformative
  periods (default 2). A crossing in force at the first lead is still detected,
  and a dip inside an otherwise degrading series no longer resets the horizon.

* **The projection mode extrapolated slopes indistinguishable from flat.**
  `min_slope_for_projection` is a magnitude gate only, so a slope of ~1e-3 with
  a one-sided p-value of 0.42 cleared it and was extrapolated. On the same
  series this produced a projected crossing 45 years out whose delta-method 95%
  interval spanned 328 years and included the past. The projection now
  additionally requires the slope to be significantly positive
  (`projection_alpha`, default 0.05); when it is not, the result is a lower
  bound, which is the honest answer for a forecast whose precision does not
  degrade in-window.

Every horizon now also carries exceedance diagnostics (`n_exceedances`,
`frac_exceedance`, `first_exceedance`, `min_run`, and `slope`/`slope_p` in trend
mode), so an isolated excursion is visible rather than silent.

New regression tests in `tests/testthat/test-shelf-life-sustained.R` pin the
defect using the real trajectory that exposed it, and assert that `min_run = 1`
and `projection_alpha = 1` reproduce the previous behaviour exactly.

## 1.3.1 — new functionality

Two new exported functions, both closing gaps identified in peer review:

* **`et_shelf_life_pool()`** pools forecast horizons across repeated forecast
  origins. Origins where the forecast was still informative when the data ran
  out are treated as right-censored rather than discarded — discarding them
  biases the pooled horizon downward, since those are exactly the long horizons.
  Estimation uses Kaplan–Meier via **survival** (a Recommended package, declared
  in `Suggests`), with a documented fallback when it is unavailable.

* **`et_priors_split()`** fits the first-stage regularized model on a subset and
  returns the disjoint complement for the Bayesian fit, so variable selection
  and the likelihood share no observations.

* **Hierarchical models** are now decomposed correctly, distinguishing
  prediction for a known group (the group effect is an estimated parameter; no
  extra variance channel) from prediction for a new level (the budget integrates
  over the group-level distribution and gains a `group_var` channel carrying the
  between-group variance).

The test suite grew from 475 to 557 tests, including numerical validation of the
decomposition for Poisson, negative-binomial and Bernoulli families against
independently simulated targets — previously only the Gaussian-identity case was
verified numerically.

## 1.3.0 — variance decomposition under autocorrelation

`decompose_uncertainty()` returned incorrect components for any model carrying a
residual autocorrelation term. It computed the environmental component as
`Var(perturbed linear predictor) - Var(unperturbed linear predictor)`. The first
array was built internally and carried no autocorrelation term; the second came
from `brms::posterior_linpred()`, which under brms's default `ar(cov = FALSE)`
parameterisation *does* include the autocorrelation contribution to the mean.
The subtraction therefore evaluated `V_env - V_AR`, which turns negative as the
AR error grows and was clamped to zero by a `pmax(., 0)` guard intended only for
Monte Carlo jitter.

The defect was invisible on inspection — the source comment asserted that "the
mean structure carries no autocorrelation" — and was found by testing against a
closed-form target. In a controlled AR(1) fit (phi = 0.9) whose analytic
environmental variance is 0.152, the reported `env_var` averaged 0.012 and was
floored to exactly zero at 13 of 15 lead times; the iid control recovered 0.173
correctly.

Fixes: both `posterior_linpred()` calls now pass `incl_autocor = FALSE` (verified
bit-identical to the internal predictor, so the AR term cancels exactly);
`var_trim` applies one variance estimator to every channel; `ar_prior` replaces
brms's flat unbounded AR prior with `normal(0, 0.5)` and warns near the unit
root; and the zero-floor is instrumented to warn when it fires more
systematically than jitter explains. Models without an autocorrelation term are
unaffected — these are provably no-ops in that case.

## Test environments

* local Linux install (Ubuntu), R 4.6.1
* win-builder (devel and release)

## R CMD check results

`R CMD check --as-cran --run-donttest` on the built tarball:

0 errors | 0 warnings | 2 notes

**Note 1 — days since last update.** This flags submission cadence, not a code
issue. If it is too soon after the previous version, I am happy to wait and
resubmit; please tell me the preferred interval.

**Note 2 — slow examples.** The expected note for a Stan-backed package:

```
* checking examples ... NOTE
Examples with CPU (user + system) or elapsed time > 5s
                        user system elapsed
decompose_uncertainty 52.230  2.827  55.159
et_calibrate          48.642  2.480  51.241
shelf_life            48.499  2.469  51.094
et_priors_split       47.540  2.533  50.203
et_fit                47.280  2.480  49.849
```

These five examples each fit one Bayesian model, so almost all of that time is
Stan model **compilation** rather than computation. The examples are already
minimal — 20 observations, one predictor, a single chain of 500 iterations, of
which sampling accounts for under a second — so the cost is not reducible by
shrinking them further. They are wrapped in `\donttest{}`.

This note was present in the previously accepted 1.2.x releases. To be precise
about what changed: the list has grown from four entries to five, because this
release adds `et_priors_split()`, whose example necessarily fits a model to show
the split workflow. No pre-existing example has been enlarged.

For context on total check time: the `testthat` suite is gated with
`skip_on_cran()`, so it runs in about 4 seconds under CRAN's settings (507 tests
locally with `NOT_CRAN=true`). Vignettes do not fit models at build time.

## Reverse dependencies

None on CRAN.
