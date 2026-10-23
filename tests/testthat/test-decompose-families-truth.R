# Non-Gaussian / non-identity-link decomposition, checked against a truth
# computed independently of the decomposition code (ErrorTracer 1.3.1).
#
# Prior to this file the decomposition was verified numerically only for the
# Gaussian-identity case (AR(1) closed form plus an iid control). Every other
# family went through the same corrected code path, which is a structural
# argument, not a numerical one. These tests supply the numerical one.
#
# Method: rather than a closed form (awkward once a non-identity link makes
# E[g^{-1}(eta)] != g^{-1}(E[eta])), the truth is obtained by BRUTE FORCE from
# the posterior draws -- simulating the response directly and taking its
# variance. That is a different computation from the analytic budget
# ErrorTracer reports, so agreement is informative rather than circular.

skip_on_cran()

# Brute-force total predictive variance: simulate the response from each draw
# and take the column variance. Uses brms::posterior_predict, which the
# decomposition does NOT use for its analytic total in the iid case.
brute_total <- function(pred) {
  apply(pred$posterior_predict, 2, stats::var)
}

# Kept as a thin alias so the per-family tests below read unchanged; the real
# check lives in expect_budget_closes() (tests/testthat/helper.R), which bounds
# the closure gap by its Monte Carlo scale now that total_var is measured rather
# than defined as the component sum.
expect_budget_reconciles <- function(dec, pred = NULL) {
  expect_budget_closes(dec, pred = pred, n_eff = 400L)
}

# The analytic total (param + residual, no env noise) should agree with the
# sampled posterior-predictive variance up to Monte Carlo error. Tolerance is
# generous because a variance estimated from n draws has relative SE
# sqrt(2/(n-1)) ~ 5% at n = 800, and the two estimators are independent.
expect_matches_brute <- function(dec, pred, tol = 0.35) {
  bt <- brute_total(pred)
  rel <- abs(dec$total_var - bt) / pmax(bt, .Machine$double.eps)
  expect_lt(stats::median(rel), tol)
}

test_that("Poisson (log link) decomposition matches a simulated truth", {
  set.seed(11)
  d <- data.frame(x = rnorm(60))
  d$y <- rpois(60, lambda = exp(1.2 + 0.6 * d$x))
  fit <- et_fit(y ~ x, data = d, family = poisson(), chains = 2, iter = 1500,
                warmup = 500, cores = 2, refresh = 0, seed = 1, silent = 2)
  nd  <- data.frame(x = seq(-1.5, 1.5, length.out = 8))

  pred <- et_predict(fit, newdata = nd, n_draws = 800, n_perturb = 400)
  dec  <- decompose_uncertainty(pred)

  expect_budget_reconciles(dec)
  expect_matches_brute(dec, pred)
  # Poisson variance function is V(mu) = mu, so residual_var must track the
  # predicted mean -- a family-specific property, not a generic one.
  mu_hat <- colMeans(pred$posterior_linpred)
  expect_gt(stats::cor(dec$residual_var, exp(mu_hat)), 0.99)
})

test_that("Negative binomial decomposition matches a simulated truth", {
  set.seed(12)
  d <- data.frame(x = rnorm(60))
  d$y <- rnbinom(60, mu = exp(1.5 + 0.5 * d$x), size = 3)
  fit <- et_fit(y ~ x, data = d, family = brms::negbinomial(), chains = 2,
                iter = 1500, warmup = 500, cores = 2, refresh = 0,
                seed = 1, silent = 2)
  nd  <- data.frame(x = seq(-1.5, 1.5, length.out = 8))

  pred <- et_predict(fit, newdata = nd, n_draws = 800, n_perturb = 400)
  dec  <- decompose_uncertainty(pred)

  expect_budget_reconciles(dec)
  expect_matches_brute(dec, pred)
  # NB residual variance mu + mu^2/phi must exceed the Poisson value mu.
  mu_hat <- exp(colMeans(pred$posterior_linpred))
  expect_true(all(dec$residual_var > mu_hat))
})

test_that("Binomial (logit link) decomposition matches a simulated truth", {
  set.seed(13)
  d <- data.frame(x = rnorm(80))
  d$y <- rbinom(80, size = 1, prob = plogis(0.3 + 1.1 * d$x))
  fit <- et_fit(y ~ x, data = d, family = brms::bernoulli(), chains = 2, iter = 1500,
                warmup = 500, cores = 2, refresh = 0, seed = 1, silent = 2)
  nd  <- data.frame(x = seq(-2, 2, length.out = 8))

  pred <- et_predict(fit, newdata = nd, n_draws = 800, n_perturb = 400)
  dec  <- decompose_uncertainty(pred)

  expect_budget_reconciles(dec)
  expect_matches_brute(dec, pred)
  # Bernoulli variance mu(1-mu) is maximised at mu = 0.5 and bounded by 0.25.
  expect_true(all(dec$residual_var <= 0.25 + 1e-8))
})

test_that("env perturbation is recovered on the response scale under a log link", {
  # With a log link the environmental channel is NOT beta^2 * delta^2 (that is
  # the link-scale answer); on the response scale it is scaled by (dmu/deta)^2 =
  # mu^2. Verify against a direct simulation rather than a link-scale formula,
  # since conflating the two scales is an easy and invisible error.
  set.seed(14)
  d <- data.frame(x = rnorm(60))
  d$y <- rpois(60, lambda = exp(1.0 + 0.8 * d$x))
  fit <- et_fit(y ~ x, data = d, family = poisson(), chains = 2, iter = 1500,
                warmup = 500, cores = 2, refresh = 0, seed = 1, silent = 2)
  nd  <- data.frame(x = seq(-1, 1, length.out = 6))
  delta <- 0.25

  dec <- decompose_uncertainty(
    et_predict(fit, newdata = nd, env_noise = list(x = delta),
               n_draws = 800, n_perturb = 800))

  expect_budget_reconciles(dec)
  expect_true(all(dec$env_var > 0))

  # Delta-method target on the response scale: Var ~= (mu * beta * delta)^2.
  dm     <- brms::as_draws_matrix(fit$fit)
  beta   <- mean(dm[, "b_x"])
  eta    <- mean(dm[, "b_Intercept"]) + beta * nd$x
  target <- (exp(eta) * beta * delta)^2
  rel    <- abs(dec$env_var - target) / target
  expect_lt(stats::median(rel), 0.5)
})

# ── 1.4.0: family-correct predictive draws and the binomial trials factor ──
# These exercise the real fitted paths that the pure-function tests in
# test-family-predictive.R cover analytically.

test_that("Gaussian identity still uses the exact location shift", {
  # The shift mu_env + (y_pp - mu) IS a draw from N(mu_env, sigma^2), so the
  # Gaussian path must be untouched by the family-RNG route. This guards the
  # case-study numbers: every analysis in the paper is Gaussian identity.
  set.seed(21)
  d <- data.frame(x = rnorm(60)); d$y <- rnorm(60, 1 + 0.8 * d$x, 0.5)
  fit <- et_fit(y ~ x, data = d, chains = 2, iter = 1000, warmup = 500,
                cores = 2, refresh = 0, seed = 1, silent = 2)
  nd <- data.frame(x = seq(-1.5, 1.5, length.out = 6))
  set.seed(3)
  p <- et_predict(fit, newdata = nd, env_noise = list(x = 0.4), n_draws = 600)

  expect_equal(
    p$predictive_draws,
    p$lp_perturbed + (p$posterior_predict - p$posterior_linpred)
  )
  expect_budget_closes(p$decomposition, pred = p)
})

test_that("Poisson predictive draws are counts and the budget closes", {
  set.seed(22)
  d <- data.frame(x = rnorm(150)); d$y <- rpois(150, exp(1.0 + 0.8 * d$x))
  fit <- et_fit(y ~ x, data = d, family = poisson(), chains = 2, iter = 1500,
                warmup = 750, cores = 2, refresh = 0, seed = 1, silent = 2)
  nd <- data.frame(x = seq(-1.5, 1.5, length.out = 8))

  set.seed(4)
  p <- et_predict(fit, newdata = nd, env_noise = list(x = 0.5), n_draws = 1500)
  dr <- p$predictive_draws

  # Pre-1.4.0 the location shift made 100% of these non-integer and 6-8%
  # negative, because it moved a Poisson(mu) residual onto mu_env.
  expect_true(all(dr >= 0))
  expect_true(all(dr == round(dr)))
  expect_budget_closes(p$decomposition, pred = p)

  # residual_var is E[V(mu_env)], so it must exceed the unperturbed E[V(mu)]
  # under a log link (the Jensen term the shift used to drop).
  set.seed(4)
  p0 <- et_predict(fit, newdata = nd, n_draws = 1500)
  expect_gt(mean(p$decomposition$residual_var),
            mean(p0$decomposition$residual_var))
})

test_that("binomial with trials > 1 decomposes on the count scale", {
  set.seed(23)
  d <- data.frame(x = rnorm(150))
  d$y <- rbinom(150, size = 10L, prob = plogis(0.3 + 1.0 * d$x))
  d$n <- 10L
  fit <- et_fit(y | trials(n) ~ x, data = d, family = binomial(), chains = 2,
                iter = 1500, warmup = 750, cores = 2, refresh = 0, seed = 1,
                silent = 2)
  nd <- data.frame(x = seq(-1.5, 1.5, length.out = 8), n = 10L)

  set.seed(5)
  p   <- et_predict(fit, newdata = nd, env_noise = list(x = 0.5), n_draws = 1500)
  dec <- p$decomposition
  dr  <- p$predictive_draws

  # Draws stay inside the support (they used to run -0.180 to 10.222).
  expect_true(all(dr >= 0 & dr <= 10))
  expect_true(all(dr == round(dr)))

  # residual_var must be n p (1 - p), not p (1 - p): with p in 0.17-0.90 and
  # n = 10 the correct values are O(1-2.5), the old ones O(0.17-0.25).
  expect_gt(mean(dec$residual_var), 1)

  # The whole budget sits on the count scale, so it closes against the measured
  # predictive variance. Pre-1.4.0 closure_gap here was ~90% of total_var.
  expect_budget_closes(dec, pred = p)

  # And driver uncertainty actually propagates: the shift left the predictive
  # variance flat in the perturbation SD.
  set.seed(5)
  p_small <- et_predict(fit, newdata = nd, env_noise = list(x = 0.1),
                        n_draws = 1500)
  expect_gt(mean(dec$total_var), 1.15 * mean(p_small$decomposition$total_var))
  expect_gt(mean(dec$env_var / dec$total_var), 0.10)
})
