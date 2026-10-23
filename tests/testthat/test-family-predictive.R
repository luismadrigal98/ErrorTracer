# tests/testthat/test-family-predictive.R
#
# Covers the two 1.4.0 family fixes (MEE-26-09-888, Reviewer 1's request to
# generalise the predictive draw beyond the Gaussian case):
#
#  A. binomial() with trials(n > 1) decomposed on the PROBABILITY scale while
#     posterior_predict() returns COUNTS. residual_var came back as p(1-p)
#     instead of n p(1-p) (n times too small), param_var as Var(p) instead of
#     Var(n p) (n^2 too small), and the env-inflated draws added a probability
#     to a count residual, landing outside the support (measured -0.180 to
#     10.222 for n = 10). The definitional total of <= 1.3.1 hid all of it;
#     the measured total reports a closure_gap of ~90% of total_var instead.
#     It survived because the only binomial test used bernoulli(size = 1),
#     where p == n p.
#
#  B. The env-inclusive predictive draws were always built as a location shift
#     mu_env + (y_pp - mu), which keeps the dispersion at mu while moving the
#     location to mu_env. Exact for Gaussian identity, wrong otherwise:
#     measured on a Poisson/log fit the predictive variance was understated by
#     0.6%, 2.4% and 4.8% at perturbation SDs 0.1 / 0.25 / 0.5, and 100% of the
#     "count" draws were non-integer with 6-8% negative. For binomial/logit the
#     variance was flat in the perturbation (2.04 at every SD) against a correct
#     2.08 -> 3.03, i.e. driver uncertainty was barely propagated at all.
#
# These tests need no Stan fit: .binomial_trials() reads the formula off a mock,
# and .family_rng() / .family_residual_var() are pure functions.

# ── A. trials resolution ────────────────────────────────────────────────────
.mock_fit <- function(family_name, lhs_expr, link = "logit") {
  list(family  = list(family = family_name, link = link),
       formula = list(formula = stats::as.formula(paste(lhs_expr, "~ x"))))
}

test_that(".binomial_trials reads trials() off the formula", {
  nd <- data.frame(x = c(-1, 0, 1), n = 10L)

  expect_equal(
    ErrorTracer:::.binomial_trials(.mock_fit("binomial", "y | trials(n)"), nd, 3L),
    c(10, 10, 10))
  # A literal is recycled to the number of rows.
  expect_equal(
    ErrorTracer:::.binomial_trials(.mock_fit("binomial", "y | trials(7)"), nd, 3L),
    c(7, 7, 7))
  # Varying trials per observation.
  nd2 <- data.frame(x = c(-1, 0, 1), n = c(2L, 5L, 9L))
  expect_equal(
    ErrorTracer:::.binomial_trials(.mock_fit("binomial", "y | trials(n)"), nd2, 3L),
    c(2, 5, 9))
})

test_that(".binomial_trials returns NULL for families that need no trials factor", {
  nd <- data.frame(x = c(-1, 0, 1), n = 10L)
  for (fam in c("gaussian", "poisson", "bernoulli", "negbinomial", "gamma")) {
    expect_null(ErrorTracer:::.binomial_trials(.mock_fit(fam, "y"), nd, 3L),
                info = fam)
  }
})

test_that(".binomial_trials errors when the trials column is absent from newdata", {
  expect_error(
    ErrorTracer:::.binomial_trials(.mock_fit("binomial", "y | trials(n)"),
                                   data.frame(x = c(-1, 0, 1)), 3L),
    regexp = "Could not resolve binomial trials"
  )
})

test_that("binomial residual_var carries the trials factor", {
  # mu is on the COUNT scale, so V(y) = mu (1 - mu / n) = n p (1 - p).
  n_s <- 50L; p <- 0.3; n_tr <- 10
  mu  <- matrix(n_tr * p, nrow = n_s, ncol = 3L)
  fam <- list(family = "binomial", link = "logit")
  disp <- list(sigma = NULL, phi = NULL, shape = NULL, nu = NULL)

  with_trials <- ErrorTracer:::.family_residual_var(fam, mu, disp,
                                                    trials = rep(n_tr, 3L))
  expect_equal(with_trials, rep(n_tr * p * (1 - p), 3L))

  # The trials factor must actually be applied. Feeding a count-scale mu with
  # trials = NULL (the combination et_predict can no longer produce) treats the
  # count as a probability and returns mu(1 - mu) = 3 * (1 - 3) < 0, which is
  # how far off the scale the two interpretations are. Up to 1.3.1 the inputs
  # were self-consistent on the PROBABILITY scale, so residual_var was p(1-p):
  # not negative, but n times too small for a count-scale total.
  without <- ErrorTracer:::.family_residual_var(fam, mu, disp, trials = NULL)
  expect_true(all(without < 0))
  expect_equal(with_trials / (n_tr * p * (1 - p)), rep(1, 3L))

  # bernoulli (n = 1) is unchanged by the fix.
  fam_b <- list(family = "bernoulli", link = "logit")
  mu_b  <- matrix(p, nrow = n_s, ncol = 3L)
  expect_equal(ErrorTracer:::.family_residual_var(fam_b, mu_b, disp),
               rep(p * (1 - p), 3L))
})

# ── B. family-correct predictive draws ──────────────────────────────────────
.disp <- function(n_s, ...) {
  d <- list(sigma = NULL, phi = NULL, shape = NULL, nu = NULL)
  utils::modifyList(d, list(...))
}

test_that(".family_rng draws inside the support and matches the variance function", {
  set.seed(5)
  n_s <- 4000L; n_obs <- 3L
  mk  <- function(v) matrix(v, nrow = n_s, ncol = n_obs, byrow = TRUE)

  # Poisson: non-negative integers, Var = mu.
  mu  <- mk(c(2, 6, 20))
  out <- ErrorTracer:::.family_rng(list(family = "poisson", link = "log"),
                                   mu, .disp(n_s))
  expect_true(all(out >= 0 & out == round(out)))
  expect_equal(apply(out, 2, var), c(2, 6, 20), tolerance = 0.12)

  # Binomial on the COUNT scale: integers in 0..n, Var = n p (1 - p).
  n_tr <- c(10, 10, 10); p <- c(0.2, 0.5, 0.9)
  out <- ErrorTracer:::.family_rng(list(family = "binomial", link = "logit"),
                                   mk(n_tr * p), .disp(n_s), trials = n_tr)
  expect_true(all(out >= 0 & out <= 10 & out == round(out)))
  expect_equal(apply(out, 2, var), n_tr * p * (1 - p), tolerance = 0.12)

  # Negative binomial: Var = mu + mu^2 / shape.
  shp <- 4
  out <- ErrorTracer:::.family_rng(
    list(family = "negbinomial", link = "log"),
    mk(c(5, 5, 5)), .disp(n_s, shape = rep(shp, n_s)))
  expect_true(all(out >= 0 & out == round(out)))
  expect_equal(mean(apply(out, 2, var)), 5 + 25 / shp, tolerance = 0.15)

  # Gaussian: the shift is exact here, but the sampler must still agree.
  out <- ErrorTracer:::.family_rng(list(family = "gaussian", link = "identity"),
                                   mk(c(0, 10, -5)), .disp(n_s, sigma = rep(2, n_s)))
  expect_equal(apply(out, 2, var), rep(4, 3L), tolerance = 0.12)
  expect_equal(colMeans(out), c(0, 10, -5), tolerance = 0.2)

  # Gamma: brms parameterises by (mu, shape), so Var = mu^2 / shape.
  out <- ErrorTracer:::.family_rng(list(family = "gamma", link = "log"),
                                   mk(c(3, 3, 3)), .disp(n_s, shape = rep(6, n_s)))
  expect_true(all(out > 0))
  expect_equal(mean(apply(out, 2, var)), 9 / 6, tolerance = 0.15)
})

test_that(".family_rng returns NULL for an unhandled family so the caller falls back", {
  mu <- matrix(0.5, nrow = 10L, ncol = 2L)
  expect_null(ErrorTracer:::.family_rng(list(family = "weibull", link = "log"),
                                        mu, .disp(10L)))
  # lognormal is deliberately absent: .family_residual_var does not support it,
  # and brms's posterior_linpred(transform = TRUE) gives the meanlog not E[y].
  expect_null(ErrorTracer:::.family_rng(list(family = "lognormal", link = "identity"),
                                        mu, .disp(10L, sigma = rep(1, 10L))))
  # Missing dispersion draws must also fall back rather than guess.
  expect_null(ErrorTracer:::.family_rng(list(family = "gaussian", link = "identity"),
                                        mu, .disp(10L)))
})

test_that("a family-correct budget closes where a location shift would not", {
  # Poisson/log fixture built by hand: mu_env is mu inflated by a log-link
  # perturbation, so E[V(mu_env)] > E[V(mu)] (the Jensen term). The budget must
  # use the mean that generated the draws, or that term is left unexplained.
  set.seed(6)
  n_s <- 6000L; n_obs <- 3L
  eta     <- matrix(rnorm(n_s * n_obs, mean = log(5), sd = 0.1), nrow = n_s)
  eps     <- matrix(rnorm(n_s * n_obs, sd = 0.4), nrow = n_s)
  mu      <- exp(eta)
  mu_env  <- exp(eta + eps)
  fam     <- list(family = "poisson", link = "log")
  disp    <- .disp(n_s)

  pp     <- matrix(rpois(n_s * n_obs, as.vector(mu)),     nrow = n_s)
  pp_env <- ErrorTracer:::.family_rng(fam, mu_env, disp)

  # Residual channel at the GENERATING mean (mu_env) -> closes.
  good <- ErrorTracer:::.decompose_from_arrays(
    pp = pp, pp_env = pp_env, mu_draws = mu, mu_perturbed = mu_env,
    mu_draws_sub = mu, family = fam, disp_draws = disp,
    mu_resid = mu_env)
  expect_lt(max(abs(good$closure_gap) / good$total_var), 0.05)

  # Residual channel at mu instead -> the Jensen term shows up as a gap.
  bad <- ErrorTracer:::.decompose_from_arrays(
    pp = pp, pp_env = pp_env, mu_draws = mu, mu_perturbed = mu_env,
    mu_draws_sub = mu, family = fam, disp_draws = disp,
    mu_resid = mu)
  expect_gt(mean(bad$closure_gap / bad$total_var),
            mean(good$closure_gap / good$total_var))
  expect_gt(mean(bad$closure_gap), 0)
})
