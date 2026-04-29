This is a resubmission. In response to the CRAN reviewer:

* Added method references with DOIs to the DESCRIPTION:
  Buerkner (2017) <doi:10.18637/jss.v080.i01> for Bayesian regression via
  Stan, Friedman, Hastie & Tibshirani (2010) <doi:10.18637/jss.v033.i01>
  for elastic net, Wright & Ziegler (2017) <doi:10.18637/jss.v077.i01>
  for random forests, and Vehtari, Gelman & Gabry (2017)
  <doi:10.1007/s11222-016-9696-4> for leave-one-out cross-validation.
* Replaced \dontrun{} with \donttest{} in all examples that need brms/Stan
  compilation (et_fit, decompose_uncertainty, shelf_life, et_calibrate);
  the wrapped examples are self-contained and execute under
  --run-donttest.
* Unwrapped the extract_priors() example (runs in <5 sec).

## Test environments
* local Linux install (Ubuntu), R 4.5.3
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 1 note

* The remaining NOTE flags this as a new submission and the
  Maintainer field — expected for an initial CRAN submission.
* `R CMD check --as-cran --run-donttest` executes all \donttest{}
  examples successfully ([322s/323s] OK).
