# Group-level variance channel (ErrorTracer 1.3.1).
#
# R1 (MEE-26-05-477, General Methods comment / Sec 4.3): "Lack of support for
# hierarchical models ... is surprising ... with hierarchical models there's an
# important distinction between in-sample and out-of-sample predictions, as the
# former uses group-specific parameters while the latter integrates over this at
# the hierarchical level."
#
# The two budgets are genuinely different objects:
#
#   in-sample  (known group): the group effect is an estimated parameter.
#                             Conditioning on it REDUCES predictive variance
#                             relative to the population level, so there is no
#                             extra channel. Var(full) - Var(pop) is NEGATIVE
#                             and must NOT be reported as a group variance.
#   out-of-sample (new group): the budget integrates over the group-level
#                             distribution, and the extra variance is tau^2.
#
# These tests pin both, and pin the closed-form target for the second.

skip_on_cran()

fit_hier <- function(tau = 1.2, sigma = 0.5, n_per = 10, n_grp = 6, seed = 1) {
  set.seed(seed)
  g  <- rep(letters[seq_len(n_grp)], each = n_per)
  u  <- stats::setNames(stats::rnorm(n_grp, 0, tau), letters[seq_len(n_grp)])
  d  <- data.frame(g = g, x = stats::rnorm(n_grp * n_per))
  d$y <- 1 + 2 * d$x + u[d$g] + stats::rnorm(nrow(d), 0, sigma)
  list(
    data = d,
    fit  = et_fit(y ~ x + (1 | g), data = d, chains = 2, iter = 1500,
                  warmup = 500, cores = 2, refresh = 0, seed = 1, silent = 2)
  )
}

test_that("a known group yields NO group channel (in-sample budget)", {
  h  <- fit_hier()
  nd <- data.frame(g = "a", x = seq(-1, 1, length.out = 5))

  pred <- et_predict(h$fit, newdata = nd, n_draws = 800, n_perturb = 400)
  dec  <- decompose_uncertainty(pred)

  # No group_var column: for a known group the effect is a parameter, and its
  # posterior uncertainty is already inside param_var.
  expect_false("group_var" %in% colnames(dec))
  # Budget still reconciles exactly.
  expect_equal(dec$total_var,
               dec$param_var + dec$env_var + dec$residual_var,
               tolerance = 1e-8)
})

test_that("a new group adds a group channel recovering tau^2", {
  h  <- fit_hier()
  nd <- data.frame(g = "brand_new_group", x = seq(-1, 1, length.out = 5))

  pred <- et_predict(h$fit, newdata = nd, n_draws = 800, n_perturb = 400)
  dec  <- decompose_uncertainty(pred)

  expect_true("group_var" %in% colnames(dec))
  expect_true(all(dec$group_var > 0))

  # The group channel should recover the fitted between-group variance. Compare
  # against tau_hat^2 from the fit itself rather than the generative tau, so the
  # test targets what the model estimated (and does not fail on sampling error
  # in tau).
  dm      <- brms::as_draws_matrix(h$fit$fit)
  tau_hat <- mean(dm[, "sd_g__Intercept"])
  expect_equal(mean(dec$group_var), tau_hat^2, tolerance = 0.5)

  # Budget reconciles exactly WITH the new channel included.
  expect_equal(dec$total_var,
               dec$param_var + dec$group_var + dec$env_var + dec$residual_var,
               tolerance = 1e-8)
})

test_that("out-of-sample prediction is wider than in-sample, as it must be", {
  h   <- fit_hier()
  x   <- seq(-1, 1, length.out = 5)
  p_in  <- et_predict(h$fit, newdata = data.frame(g = "a", x = x),
                      n_draws = 800, n_perturb = 400)
  p_out <- et_predict(h$fit, newdata = data.frame(g = "new_grp", x = x),
                      n_draws = 800, n_perturb = 400)

  d_in  <- decompose_uncertainty(p_in)
  d_out <- decompose_uncertainty(p_out)

  # Integrating over the group level must cost variance. This is the substantive
  # statement: reporting the in-sample budget for a new group would understate
  # forecast uncertainty by roughly tau^2.
  expect_true(mean(d_out$total_var) > mean(d_in$total_var))
})

test_that("a model without group terms is completely unaffected", {
  set.seed(2)
  d <- data.frame(x = rnorm(40))
  d$y <- 1 + 2 * d$x + rnorm(40, 0, 0.5)
  fit <- et_fit(y ~ x, data = d, chains = 2, iter = 1000, warmup = 500,
                cores = 2, refresh = 0, seed = 1, silent = 2)
  nd  <- data.frame(x = seq(-1, 1, length.out = 5))

  dec <- decompose_uncertainty(et_predict(fit, newdata = nd, n_draws = 600,
                                          n_perturb = 300))
  expect_false("group_var" %in% colnames(dec))
  expect_equal(dec$total_var,
               dec$param_var + dec$env_var + dec$residual_var,
               tolerance = 1e-8)
})
