This is a minor release (1.0.2 -> 1.1.0) adding new functionality and
one bug fix. No reverse dependencies were detected on CRAN.

## Summary of changes

* New: response-scale variance decomposition for non-Gaussian families
  (Binomial, Poisson, Student-t, Negative-Binomial, Beta, Gamma) via the
  inverse link. Gaussian identity is numerically unchanged.
* New: `et_predict()` gains `n_env_draws` for averaging multiple
  environmental perturbations per posterior draw; the decomposition table
  gains a `v_env_mcse` column reporting the chi-squared standard error of
  `env_var`.
* New: `decompose_uncertainty()` reports a fourth `temporal_var`
  component when the model formula contains an autocorrelation term
  (`ar()`, `ma()`, `arma()`, `cosy()`, `unstr()`, `sar()`, `car()`).
  `et_plot_decomposition()` and `print.et_prediction()` surface the new
  component automatically.
* New: `shelf_life()` now returns `se_t_star` in the projected-mode
  horizon attribute (delta-method SE of `t* = (tau - a) / b`).
* Deprecation: `shelf_life(plausible_range = ...)` is deprecated in
  favour of `response_scale`. The old argument still works and emits a
  warning; it will be removed in a future release.
* Fix: `et_plot_calibration()` now auto-detects the grouping column when
  it is named something other than `group` (e.g. `cluster_id`, `species`,
  `regime`). Adds a `group_col` argument for explicit control, including
  `group_col = NA` to force a single un-grouped series.

## Test environments

* local Linux install (Ubuntu), R 4.5.3
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

`R CMD check --as-cran --run-donttest` executes all \donttest{}
examples successfully.

## Reverse dependencies

None on CRAN (checked with `revdepcheck::revdep_check()`).
