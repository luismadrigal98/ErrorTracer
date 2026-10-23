# tests/testthat/test-draw-pairing.R
#
# Regression tests for the draw-pairing fix and the measured-total
# decomposition introduced in 1.4.0 (MEE-26-09-888, Reviewer 1).
#
# The defect: .inflate_env_predictive() used to RECYCLE a shorter mu_perturbed
# across all n_draws rows,
#     idx <- rep_len(seq_len(n_perturb), n_draws)
#     mu_perturbed[idx, ] + (pp - mu_draws)
# which paired the env-perturbed mean of draw `s mod n_perturb` with the
# residual of draw `s`. That destroys Cov_s(mu^(s), y^pp,(s) - mu^(s)). The
# covariance is ~0 for an iid fit but substantially NEGATIVE for an ar() fit,
# where a draw with large |phi| carries a wide accumulated residual, so the
# predictive variance was inflated. Measured on the Portal rodent series the
# gap V(predictive) - V(pp) was 20x, 22x and 196x env_var for the three most
# persistent species, against 0.75-1.33x under matched pairing. Those draws are
# the single source of truth for the credible interval, PIT, calibration AND
# shelf life, so the inflation reached all four.
#
# None of these tests needs a Stan fit: the fixture below reproduces the
# negative mean/residual covariance analytically.

# ── Fixture: known variances and a strong negative Cov(mu, residual) ─────────
# z is a per-draw common factor; the residual loads on it with the opposite
# sign, which is what an AR posterior does (large |phi| => wide residual).
#
#   Var(mu)            = 1.00     -> param_var
#   Var(resid)         = 1.25     -> residual_var (sigma fixed at sqrt(1.25))
#   Cov(mu, resid)     = -1.00
#   Var(pp)            = 1 + 1.25 - 2 = 0.25
#   Var(perturbation)  = 0.09     -> env_var
#   Var(pp_env)        = 0.25 + 0.09 = 0.34
.pair_fixture <- function(n_draws = 8000L, n_obs = 3L, seed = 7L) {
  set.seed(seed)
  z <- rnorm(n_draws)
  w <- matrix(rnorm(n_draws * n_obs), nrow = n_draws)
  p <- matrix(rnorm(n_draws * n_obs, sd = 0.3), nrow = n_draws)
  means <- matrix(rep(seq_len(n_obs), each = n_draws), nrow = n_draws)

  mu    <- means + z                      # Var 1, shared factor
  resid <- -z + 0.5 * w                   # Var 1.25, Cov(mu, resid) = -1
  list(
    pp           = mu + resid,
    mu_draws     = mu,
    mu_perturbed = mu + p,                # env_var = Var(p) = 0.09
    sigma        = rep(sqrt(1.25), n_draws),
    n_obs        = n_obs
  )
}

.gauss_family <- list(family = "gaussian", link = "identity",
                      linkinv = identity)

.decomp_fixture <- function(f, has_autocor = TRUE) {
  pp_env <- ErrorTracer:::.inflate_env_predictive(f$mu_perturbed, f$pp,
                                                  f$mu_draws)
  ErrorTracer:::.decompose_from_arrays(
    pp           = f$pp,
    pp_env       = pp_env,
    mu_draws     = f$mu_draws,
    mu_perturbed = f$mu_perturbed,
    mu_draws_sub = f$mu_draws,
    family       = .gauss_family,
    disp_draws   = list(sigma = f$sigma, phi = NULL, shape = NULL, nu = NULL),
    has_autocor  = has_autocor
  )
}

# ── 1. The helper refuses to recycle ────────────────────────────────────────
test_that(".inflate_env_predictive requires matched draw counts", {
  f <- .pair_fixture(n_draws = 400L)
  half <- f$mu_perturbed[seq_len(200L), , drop = FALSE]
  expect_error(
    ErrorTracer:::.inflate_env_predictive(half, f$pp, f$mu_draws),
    regexp = "nrow"
  )
  expect_silent(
    ErrorTracer:::.inflate_env_predictive(f$mu_perturbed, f$pp, f$mu_draws)
  )
})

# ── 2. Matched pairing recovers env_var; recycling does not ─────────────────
test_that("matched pairing recovers env_var, recycling inflates it", {
  f <- .pair_fixture()
  n_draws <- nrow(f$pp)

  v_pp     <- apply(f$pp, 2, var)
  env_true <- apply(f$mu_perturbed, 2, var) - apply(f$mu_draws, 2, var)

  # Matched: V(pp_env) - V(pp) is the environmental increment. Both sides are
  # differences of two nearly equal sample variances, so each column carries
  # O(var * sqrt(2 / n_draws)) of Monte Carlo error -- compare the means, and
  # let the ratio assertions below carry the discrimination.
  matched <- ErrorTracer:::.inflate_env_predictive(f$mu_perturbed, f$pp,
                                                   f$mu_draws)
  gap_matched <- apply(matched, 2, var) - v_pp
  expect_equal(mean(gap_matched), mean(env_true), tolerance = 0.15)

  # Recycled (the pre-1.4.0 formula, reproduced here so the test fails if it
  # ever comes back): a quarter of the rows keep their own residual, so the
  # destroyed covariance shows up as a large positive bias.
  idx          <- rep_len(seq_len(n_draws / 4L), n_draws)
  recycled     <- f$mu_perturbed[idx, , drop = FALSE] + (f$pp - f$mu_draws)
  gap_recycled <- apply(recycled, 2, var) - v_pp

  expect_gt(mean(gap_recycled / env_true), 10)
  expect_lt(mean(gap_matched  / env_true), 1.2)
})

# ── 3. total_var is measured, and the channels close on it ─────────────────
test_that("total_var equals the env-inclusive predictive variance", {
  f <- .pair_fixture()
  d <- .decomp_fixture(f, has_autocor = TRUE)

  pp_env <- ErrorTracer:::.inflate_env_predictive(f$mu_perturbed, f$pp,
                                                  f$mu_draws)
  expect_equal(d$total_var, apply(pp_env, 2, var))

  # For an autocorrelation model temporal_var is the remainder, so the budget
  # closes exactly rather than by definition.
  expect_equal(
    d$param_var + d$env_var + d$residual_var + d$temporal_var,
    d$total_var
  )
  expect_false("closure_gap" %in% colnames(d))
})

test_that("iid models report a small closure_gap instead of a temporal channel", {
  f <- .pair_fixture()
  d <- .decomp_fixture(f, has_autocor = FALSE)

  expect_false("temporal_var" %in% colnames(d))
  expect_true("closure_gap" %in% colnames(d))
  expect_equal(
    d$closure_gap,
    d$total_var - (d$param_var + d$env_var + d$residual_var)
  )
})

# ── 4. temporal_var is signed, not clamped at zero ─────────────────────────
test_that("temporal_var stays negative when errors reduce predictive variance", {
  f <- .pair_fixture()
  d <- .decomp_fixture(f, has_autocor = TRUE)

  # Cov(mu, resid) = -1 by construction, so the predictive variance (0.34) is
  # far BELOW the independent-innovation baseline (1 + 0.09 + 1.25 = 2.34).
  # Pre-1.4.0 pmax(., 0) reported 0 here and silently broke the budget.
  expect_lt(mean(d$temporal_var), -1.5)
  expect_equal(mean(d$param_var),    1.00, tolerance = 0.05)
  expect_equal(mean(d$residual_var), 1.25, tolerance = 1e-8)
  expect_equal(mean(d$env_var),      0.09, tolerance = 0.05)
})

# ── 5. Zero env noise leaves the total equal to V(pp) ─────────────────────
test_that("with no perturbation total_var reduces to the plain predictive variance", {
  f <- .pair_fixture()
  f$mu_perturbed <- f$mu_draws          # zero-noise path
  d <- .decomp_fixture(f, has_autocor = TRUE)

  expect_equal(d$total_var, apply(f$pp, 2, var))
  expect_equal(d$env_var, rep(0, f$n_obs))
  expect_equal(d$v_env_mcse, rep(0, f$n_obs))
})
