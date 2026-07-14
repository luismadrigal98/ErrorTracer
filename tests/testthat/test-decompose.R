# tests/testthat/test-decompose.R

# Minimal mock et_prediction (Gaussian identity) for tests that don't need
# a real Stan fit.  All decomposition is computed via the internal helpers so
# that the mock stays in sync with the production code path.
.mock_et_prediction <- function(n_obs = 5, n_draws = 100) {
  set.seed(99)
  pp <- matrix(rnorm(n_draws * n_obs, mean = 0.5), nrow = n_draws)
  lp <- matrix(rnorm(n_draws * n_obs, mean = 0.5, sd = 0.3), nrow = n_draws)
  # Perturbed LP — slightly higher variance than lp; for Gaussian identity
  # mu_draws == lp and mu_perturbed == lp_p.
  lp_p <- lp + matrix(rnorm(n_draws * n_obs, sd = 0.1), nrow = n_draws)
  sigma_draws <- abs(rnorm(n_draws, mean = 0.5, sd = 0.1))

  gauss_family <- list(family = "gaussian", link = "identity",
                       linkinv = identity)
  disp_draws   <- list(sigma = sigma_draws, phi = NULL,
                       shape = NULL, nu = NULL)

  decomp <- ErrorTracer:::.decompose_from_arrays(
    pp           = pp,
    mu_draws     = lp,       # Gaussian: mu_draws == lp
    mu_perturbed = lp_p,
    mu_draws_sub = lp,       # n_perturb == n_draws in mock
    family       = gauss_family,
    disp_draws   = disp_draws
  )

  ci_df <- ErrorTracer:::.compute_ci(pp, c(0.5, 0.9, 0.95))

  structure(
    list(
      posterior_predict  = pp,
      posterior_linpred  = lp,
      lp_perturbed       = lp_p,
      sigma_draws        = sigma_draws,
      credible_intervals = ci_df,
      decomposition      = decomp,
      newdata            = data.frame(obs_id = seq_len(n_obs)),
      model              = NULL,
      env_noise          = NULL,
      n_draws            = n_draws
    ),
    class = "et_prediction"
  )
}

test_that("decompose_uncertainty returns correct columns", {
  pred   <- .mock_et_prediction()
  decomp <- decompose_uncertainty(pred)

  expect_s3_class(decomp, "data.frame")
  expect_true(all(c("obs_id", "param_var", "env_var", "v_env_mcse",
                    "residual_var", "total_var") %in% colnames(decomp)))
})

test_that("decompose_uncertainty has non-negative variance components", {
  pred   <- .mock_et_prediction(n_obs = 10)
  decomp <- decompose_uncertainty(pred)

  expect_true(all(decomp$param_var    >= 0))
  expect_true(all(decomp$env_var      >= 0))
  expect_true(all(decomp$v_env_mcse   >= 0))
  expect_true(all(decomp$residual_var >= 0))
  expect_true(all(decomp$total_var    >= 0))
})

test_that("decompose_uncertainty row count matches n_obs", {
  n_obs <- 7
  pred  <- .mock_et_prediction(n_obs = n_obs)
  decomp <- decompose_uncertainty(pred)
  expect_equal(nrow(decomp), n_obs)
})

test_that("decompose_uncertainty.default raises an error", {
  expect_error(decompose_uncertainty(list(a = 1)), "et_prediction")
})

test_that("param_var < total_var on average (residual inflates total)", {
  pred   <- .mock_et_prediction(n_obs = 20, n_draws = 500)
  decomp <- decompose_uncertainty(pred)
  expect_true(mean(decomp$total_var) > mean(decomp$param_var))
})

test_that("residual_var is constant across observations for Gaussian", {
  # For Gaussian identity sigma is constant across obs (scalar sigma^2)
  pred   <- .mock_et_prediction(n_obs = 6)
  decomp <- decompose_uncertainty(pred)
  expect_equal(length(unique(round(decomp$residual_var, 10))), 1L)
})

test_that("v_env_mcse is finite and positive when env noise is present", {
  pred   <- .mock_et_prediction(n_obs = 5, n_draws = 50)
  decomp <- decompose_uncertainty(pred)
  # Some env noise was injected in the mock (lp_p = lp + noise)
  expect_true(all(is.finite(decomp$v_env_mcse)))
  expect_true(all(decomp$v_env_mcse >= 0))
})

test_that("temporal_var column is absent unless has_autocor = TRUE", {
  pred   <- .mock_et_prediction()
  decomp <- decompose_uncertainty(pred)
  expect_false("temporal_var" %in% colnames(decomp))
})

test_that("temporal_var captures the autocor gap when has_autocor = TRUE", {
  # Build a mock where pp has variance much greater than param + env + residual,
  # mimicking AR-induced predictive variance growth.
  set.seed(101)
  n_obs   <- 5
  n_draws <- 200
  # Tight linpred (small param), perturbed lp ~ lp (small env), residual sigma ~ 0.3,
  # but pp inflated by extra AR-style innovation (sd = 0.8) on top of lp.
  lp   <- matrix(rnorm(n_draws * n_obs, mean = 0.5, sd = 0.05), nrow = n_draws)
  lp_p <- lp + matrix(rnorm(n_draws * n_obs, sd = 0.01), nrow = n_draws)
  sigma_draws <- rep(0.3, n_draws)
  pp <- lp + matrix(rnorm(n_draws * n_obs, sd = sqrt(0.3^2 + 0.8^2)), nrow = n_draws)

  gauss_family <- list(family = "gaussian", link = "identity",
                       linkinv = identity)
  disp_draws   <- list(sigma = sigma_draws, phi = NULL,
                       shape = NULL, nu = NULL)

  d_iid <- ErrorTracer:::.decompose_from_arrays(
    pp = pp, mu_draws = lp, mu_perturbed = lp_p, mu_draws_sub = lp,
    family = gauss_family, disp_draws = disp_draws, has_autocor = FALSE)
  d_ar <- ErrorTracer:::.decompose_from_arrays(
    pp = pp, mu_draws = lp, mu_perturbed = lp_p, mu_draws_sub = lp,
    family = gauss_family, disp_draws = disp_draws, has_autocor = TRUE)

  expect_false("temporal_var" %in% colnames(d_iid))
  expect_true("temporal_var" %in% colnames(d_ar))
  expect_true(all(d_ar$temporal_var >= 0))
  # total_var is DEFINED as the sum of every component, so param + env +
  # residual + temporal reconstructs total_var EXACTLY (env now enters the
  # total instead of being bolted on top).
  recon <- d_ar$param_var + d_ar$env_var + d_ar$residual_var + d_ar$temporal_var
  expect_equal(recon, d_ar$total_var, tolerance = 1e-10)
  # And the temporal gap should be substantial relative to residual_var alone.
  expect_true(mean(d_ar$temporal_var) > mean(d_ar$residual_var))
})

test_that("autocorrelation temporal_var is robust to a single explosive draw", {
  # Regression guard for the near-unit-root pathology: when an AR posterior has
  # a tiny fraction of |phi| >= 1 draws, the far-lead forecast variance is
  # heavy-tailed and a raw var() would be dominated by one explosive draw while
  # the quantile interval is unaffected. The decomposition now uses a winsorized
  # (tail-robust) predictive variance so it stays consistent with that interval.
  set.seed(202)
  n_obs <- 4L; n_draws <- 500L
  lp    <- matrix(rnorm(n_draws * n_obs, mean = 0.5, sd = 0.05), nrow = n_draws)
  sigma_draws <- rep(0.3, n_draws)
  # Well-behaved AR-style pp: linpred + a moderate innovation.
  pp <- lp + matrix(rnorm(n_draws * n_obs, sd = 0.4), nrow = n_draws)
  # Inject ONE explosive (near-non-stationary) draw in the last column.
  pp_bad <- pp
  pp_bad[1L, n_obs] <- pp_bad[1L, n_obs] + 60

  gauss_family <- list(family = "gaussian", link = "identity",
                       linkinv = identity)
  disp_draws   <- list(sigma = sigma_draws, phi = NULL, shape = NULL, nu = NULL)

  d_clean <- ErrorTracer:::.decompose_from_arrays(
    pp = pp, mu_draws = lp, mu_perturbed = lp, mu_draws_sub = lp,
    family = gauss_family, disp_draws = disp_draws, has_autocor = TRUE)
  d_bad <- suppressWarnings(ErrorTracer:::.decompose_from_arrays(
    pp = pp_bad, mu_draws = lp, mu_perturbed = lp, mu_draws_sub = lp,
    family = gauss_family, disp_draws = disp_draws, has_autocor = TRUE))

  # The one explosive draw must not blow up temporal_var: the robust estimate
  # stays close to the clean one (a raw var() would be several times larger).
  expect_lt(d_bad$temporal_var[n_obs], 2 * d_clean$temporal_var[n_obs])
  # Budget still reconciles exactly.
  recon <- d_bad$param_var + d_bad$env_var + d_bad$residual_var + d_bad$temporal_var
  expect_equal(recon, d_bad$total_var, tolerance = 1e-10)
})

test_that("iid budget reconciles exactly: param + env + residual == total", {
  # T1: total_var must CONTAIN env_var so the percentage shares sum to 100%.
  set.seed(202)
  n_obs <- 8; n_draws <- 400
  lp   <- matrix(rnorm(n_draws * n_obs, mean = 0.5, sd = 0.25), nrow = n_draws)
  # Perturbed mean has clearly larger spread → non-trivial env_var.
  lp_p <- lp + matrix(rnorm(n_draws * n_obs, sd = 0.4), nrow = n_draws)
  sigma_draws <- abs(rnorm(n_draws, mean = 0.6, sd = 0.05))
  pp   <- lp + matrix(rnorm(n_draws * n_obs, sd = rep(sigma_draws, n_obs)),
                      nrow = n_draws)

  gauss_family <- list(family = "gaussian", link = "identity",
                       linkinv = identity)
  disp_draws   <- list(sigma = sigma_draws, phi = NULL, shape = NULL, nu = NULL)

  d <- ErrorTracer:::.decompose_from_arrays(
    pp = pp, mu_draws = lp, mu_perturbed = lp_p, mu_draws_sub = lp,
    family = gauss_family, disp_draws = disp_draws, has_autocor = FALSE)

  # env_var is genuinely non-zero here
  expect_true(mean(d$env_var) > 0)
  # Exact reconciliation → percentage shares sum to 100%
  recon <- d$param_var + d$env_var + d$residual_var
  expect_equal(recon, d$total_var, tolerance = 1e-10)
  pct <- 100 * recon / d$total_var
  expect_equal(pct, rep(100, n_obs), tolerance = 1e-8)
})

test_that(".formula_has_autocor detects ar()/ma()/arma() and rejects plain formulas", {
  expect_false(ErrorTracer:::.formula_has_autocor(y ~ x))
  expect_true(ErrorTracer:::.formula_has_autocor(y ~ x + ar(time = t, p = 1)))
  expect_true(ErrorTracer:::.formula_has_autocor(y ~ x + ma(time = t, q = 1)))
})

test_that("decompose_uncertainty.et_prediction_list returns grouped data.frame", {
  p1 <- .mock_et_prediction(n_obs = 4)
  p2 <- .mock_et_prediction(n_obs = 3)

  pred_list <- structure(
    list(
      predictions = list(A = p1, B = p2),
      grouping    = "cluster_id",
      newdata     = data.frame(cluster_id = rep(c("A","B"), c(4,3)))
    ),
    class = "et_prediction_list"
  )

  decomp <- decompose_uncertainty(pred_list)
  expect_true("group" %in% colnames(decomp))
  expect_equal(nrow(decomp), 7L)
  expect_equal(sort(unique(decomp$group)), c("A", "B"))
})

test_that("decompose_uncertainty.et_prediction_list skips NULL groups", {
  p1 <- .mock_et_prediction(n_obs = 4)

  pred_list <- structure(
    list(
      predictions = list(A = p1, B = NULL),
      grouping    = "cluster_id",
      newdata     = data.frame(cluster_id = rep("A", 4))
    ),
    class = "et_prediction_list"
  )

  decomp <- decompose_uncertainty(pred_list)
  expect_equal(nrow(decomp), 4L)
  expect_false("B" %in% decomp$group)
})
