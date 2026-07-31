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

expect_budget_reconciles <- function(dec) {
  comp <- dec$param_var + dec$env_var + dec$residual_var
  if ("group_var" %in% colnames(dec))    comp <- comp + dec$group_var
  if ("temporal_var" %in% colnames(dec)) comp <- comp + dec$temporal_var
  expect_equal(dec$total_var, comp, tolerance = 1e-8)
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
