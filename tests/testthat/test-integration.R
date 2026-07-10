# tests/testthat/test-integration.R
#
# Integration tests that exercise et_fit(), et_predict(), et_calibrate(), and
# et_diagnose() with a real (but minimal) Stan/brms fit.
#
# These tests are skipped on CRAN (Stan compilation is too slow) and require
# the rstan backend.  They are the primary source of coverage for:
#   R/fit.R, R/calibrate.R, and the brms-calling branches of R/predict.R
#
# Runtime: ~60-120 s on first call (C++ compilation); <10 s on reruns.

# ── Shared fixture: fit a minimal model once per file ───────────────────────

# testthat caches the environment per file, so we compute the fit at top level.
skip_if_not_installed("brms")
# covr 3.6.5 does not set NOT_CRAN, so skip_on_cran() would skip these tests
# during coverage measurement.  Detect the covr subprocess and opt in.
if (requireNamespace("covr", quietly = TRUE) && covr::in_covr()) {
  Sys.setenv(NOT_CRAN = "true")
}
skip_on_cran()

options(brms.backend = "rstan")

# Use a tiny synthetic dataset (n = 15) for speed
set.seed(42L)
.n   <- 15L
.int_df <- data.frame(
  y  = rnorm(.n, mean = 0.5),
  x1 = rnorm(.n),
  x2 = rnorm(.n)
)
.int_valid <- data.frame(
  y  = rnorm(5L, mean = 0.5),
  x1 = rnorm(5L),
  x2 = rnorm(5L)
)

# Extract priors via lm (no glmnet needed)
.lm_fit    <- lm(y ~ x1 + x2, data = .int_df)
.int_prior <- extract_priors(.lm_fit, multiplier = 2.0, min_sd = 0.1)

# Fit once — used by all tests in this file
suppressWarnings(
  .int_fit <- et_fit(
    formula = y ~ x1 + x2,
    data    = .int_df,
    priors  = .int_prior,
    chains  = 1L, iter = 400L, warmup = 200L,
    cores   = 1L, seed  = 42L, refresh = 0L,
    silent  = 2L
  )
)

# ── et_fit tests ─────────────────────────────────────────────────────────────

test_that("et_fit returns an et_model with the correct structure", {
  expect_s3_class(.int_fit, "et_model")
  expect_true(inherits(.int_fit$fit, "brmsfit"))
  expect_equal(nrow(.int_fit$data), .n)
  expect_true(!is.null(.int_fit$prior_spec))
})

test_that("print.et_model runs without error", {
  expect_output(print(.int_fit), "et_model")
})

test_that("summary.et_model runs without error", {
  expect_output(summary(.int_fit), "Fixed effects")
})

test_that("et_fit Rhat is below 1.10 for trivial model", {
  rhats <- brms::rhat(.int_fit$fit)
  expect_true(all(rhats < 1.10, na.rm = TRUE))
})

# ── et_predict tests ─────────────────────────────────────────────────────────

test_that("et_predict returns an et_prediction with correct dimensions", {
  pred <- et_predict(
    model     = .int_fit,
    newdata   = .int_valid,
    env_noise = list(x1 = 0.2, x2 = 0.1),
    n_draws   = 200L,
    ci_levels = c(0.50, 0.90),
    n_perturb = 50L
  )
  expect_s3_class(pred, "et_prediction")
  expect_equal(ncol(pred$posterior_predict), 5L)
  expect_equal(nrow(pred$posterior_predict), 200L)
  expect_true(all(c("row_id", "ci_level", "lower", "upper", "width") %in%
                    colnames(pred$credible_intervals)))
})

test_that("et_predict with no env_noise gives zero env_var", {
  pred_no_noise <- suppressWarnings(et_predict(
    model     = .int_fit,
    newdata   = .int_valid,
    env_noise = NULL,
    n_draws   = 100L,
    n_perturb = 50L
  ))
  expect_true(all(pred_no_noise$decomposition$env_var == 0))
})

test_that("et_predict with time-varying env_noise inflates later env_var", {
  noise_vec <- c(0.01, 0.01, 0.5, 0.5, 0.5)   # low for obs 1-2, high for 3-5
  pred_tv <- suppressWarnings(et_predict(
    model     = .int_fit,
    newdata   = .int_valid,
    env_noise = list(x1 = noise_vec),
    n_draws   = 200L,
    n_perturb = 200L
  ))
  env_var <- pred_tv$decomposition$env_var
  # High-noise obs (3-5) should have larger env_var than low-noise (1-2)
  expect_gt(mean(env_var[3:5]), mean(env_var[1:2]))
})

# ── et_diagnose tests ────────────────────────────────────────────────────────

test_that("et_diagnose returns convergence summary without LOO", {
  diag <- et_diagnose(.int_fit, loo = FALSE)
  expect_true(!is.null(diag$convergence))
  expect_true(is.logical(diag$convergence$rhat_all_ok))
  expect_true(is.logical(diag$convergence$neff_all_ok))
  expect_true(is.numeric(diag$convergence$n_divergences))
  expect_null(diag$loo)
})

test_that("et_diagnose with loo = TRUE returns LOO summary", {
  diag <- suppressWarnings(et_diagnose(.int_fit, loo = TRUE))
  expect_false(is.null(diag$loo))
  expect_true(is.numeric(diag$loo$elpd_loo))
})

# ── et_calibrate tests ───────────────────────────────────────────────────────

test_that("et_calibrate returns coverage data.frame with correct structure", {
  pred <- suppressWarnings(et_predict(
    .int_fit, .int_valid, n_draws = 200L, n_perturb = 50L
  ))
  cal <- et_calibrate(pred, observed = .int_valid,
                      response_col = "y", ci_levels = c(0.50, 0.90))
  expect_s3_class(cal, "data.frame")
  expect_true(all(c("ci_level", "nominal", "observed_coverage",
                    "n_obs", "calibration_error") %in% colnames(cal)))
  expect_equal(nrow(cal), 2L)
  expect_true(all(cal$observed_coverage >= 0 & cal$observed_coverage <= 1))
})

test_that("et_calibrate print method runs without error", {
  pred <- suppressWarnings(et_predict(
    .int_fit, .int_valid, n_draws = 200L, n_perturb = 50L
  ))
  cal <- et_calibrate(pred, observed = .int_valid, response_col = "y")
  expect_output(print(cal), "calibration")
})

# ── AR forecast accumulation (T3) ────────────────────────────────────────────

test_that("AR(1) forecast variance accumulates with lead time", {
  # brms accumulates AR variance across a horizon only when the forecast rows
  # continue the training series with response = NA (see
  # .posterior_predict_forecast). This guards against the flat-CI regression
  # the reviewer flagged (L335).
  set.seed(7L)
  nn <- 50L; ph <- 0.85
  ee <- numeric(nn); ee[1] <- rnorm(1)
  for (tt in 2:nn) ee[tt] <- ph * ee[tt - 1] + rnorm(1)
  xx <- rnorm(nn)
  ar_df <- data.frame(t = seq_len(nn), x = xx, y = 1 + 0.5 * xx + ee)
  tr <- ar_df[1:35, ]; fu <- ar_df[36:50, ]

  ar_fit <- suppressWarnings(et_fit(
    y ~ x + ar(time = t, p = 1), data = tr,
    chains = 1L, iter = 600L, warmup = 300L, cores = 1L, seed = 7L, refresh = 0L
  ))
  pr <- suppressWarnings(et_predict(ar_fit, newdata = fu, n_draws = 300L))
  d  <- decompose_uncertainty(pr)
  H  <- nrow(fu)

  expect_true("temporal_var" %in% colnames(d))
  # Total predictive variance grows substantially from lead 1 to lead H
  expect_gt(d$total_var[H], 1.5 * d$total_var[1])
  # So does the reported 95% CI width
  w95 <- pr$credible_intervals
  w95 <- w95$width[w95$ci_level == 0.95]
  expect_gt(w95[H], 1.3 * w95[1])
  # The 4-way budget still reconciles exactly
  recon <- d$param_var + d$env_var + d$residual_var + d$temporal_var
  expect_equal(recon, d$total_var, tolerance = 1e-8)
})

# ── et_fit with lm prior (default method path) ───────────────────────────────

test_that("et_fit works with flat (NULL) priors", {
  suppressWarnings(
    fit_flat <- et_fit(
      formula = y ~ x1 + x2,
      data    = .int_df,
      priors  = NULL,
      chains  = 1L, iter = 300L, warmup = 150L,
      cores   = 1L, seed = 7L, refresh = 0L, silent = 2L
    )
  )
  expect_s3_class(fit_flat, "et_model")
})
