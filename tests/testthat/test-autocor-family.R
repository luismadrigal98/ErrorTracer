# tests/testthat/test-autocor-family.R
#
# Regression tests for the non-Gaussian autocorrelation gate added in 1.4.3.
#
# The defect: brms represents an ar() term two different ways. For gaussian and
# student it uses "natural" residuals (err[n] = Y[n] - mu[n]) and can iterate
# the process onto new rows. For every other family with the DEFAULT
# cov = FALSE it declares a latent parameter vector instead
# (err = sderr * zerr, added to the linear predictor) and then refuses to
# predict: "Cannot predict new latent residuals when using cov = FALSE in
# autocor terms."
#
# Before the gate, et_fit(y ~ x + ar(time = t), family = poisson()) compiled and
# sampled happily and et_predict() died several frames inside brms with that
# message and no remedy. Measured on a 45-observation Poisson fit.
#
# With cov = TRUE prediction works, but the env-inclusive draws keep the
# location-shift construction (see .inflate_env_predictive) and are therefore
# not guaranteed to respect the response support. Measured on the same data:
# 12.3% of draws negative, 100% non-integer, 41 of 60 credible-interval rows
# with a negative lower bound. That is a warning, not an error, because the
# construction is deliberate -- y_pp - mu is what carries the forward iteration.
#
# Most of this file needs no Stan: the gate is formula-and-family algebra.

# ── .family_has_natural_residuals(): gate on the family NAME, not the link ──
# brms computes err = Y - mu for gaussian and student whatever the link, so a
# log-link gaussian is still predictable.
test_that(".family_has_natural_residuals() follows brms's residual handling", {
  expect_true(.family_has_natural_residuals(NULL))              # gaussian default
  expect_true(.family_has_natural_residuals(stats::gaussian()))
  expect_true(.family_has_natural_residuals(stats::gaussian(link = "log")))
  expect_true(.family_has_natural_residuals(brms::student()))
  expect_false(.family_has_natural_residuals(stats::poisson()))
  expect_false(.family_has_natural_residuals(stats::binomial()))
  expect_false(.family_has_natural_residuals(brms::negbinomial()))
})

# ── .autocor_cov(): read cov off the term via brms's own constructors ───────
test_that(".autocor_cov() reports the brms cov parameterisation", {
  expect_false(.autocor_cov(brms::bf(y ~ x + ar(time = t))))
  expect_true(.autocor_cov(brms::bf(y ~ x + ar(time = t, cov = TRUE))))
  # brms rejects ar(p = 2, cov = TRUE) outright ("Covariance formulation of
  # ARMA structures is only possible for effects of maximal order one"), so the
  # term cannot be built and NA is the honest answer. The practical
  # consequence, which the error message states: a non-Gaussian family cannot
  # carry an autocorrelation term beyond lag 1 in brms at all.
  expect_true(is.na(.autocor_cov(brms::bf(y ~ x + ar(time = t, p = 2, cov = TRUE)))))
  expect_false(.autocor_cov(brms::bf(y ~ x + ar(time = t, p = 2))))
  # No autocorrelation term at all: nothing to report.
  expect_true(is.na(.autocor_cov(brms::bf(y ~ x))))
  expect_true(is.na(.autocor_cov(NULL)))
})

test_that(".autocor_cov() handles ma() and arma() terms", {
  # ma()/arma() require cov = TRUE in brms regardless of family, so the reader
  # must still return the right answer for them.
  expect_true(.autocor_cov(brms::bf(y ~ x + ma(time = t, cov = TRUE))))
  expect_true(.autocor_cov(brms::bf(y ~ x + arma(time = t, cov = TRUE))))
})

# ── .check_autocor_family(): errors only where brms cannot predict ──────────
test_that(".check_autocor_family() blocks non-Gaussian ar(cov = FALSE)", {
  f <- brms::bf(y ~ x + ar(time = t))
  expect_error(.check_autocor_family(f, stats::poisson(), "fit"),
               "requires cov = TRUE")
  expect_error(.check_autocor_family(f, stats::poisson(), "fit"),
               "latent residuals")
  # The message has to carry the remedy, not just the diagnosis.
  expect_error(.check_autocor_family(f, stats::poisson(), "fit"),
               "cov = TRUE\\)")
  expect_error(.check_autocor_family(f, stats::poisson(), "fit"),
               "order one only")
  # Wording differs by stage but the gate is the same.
  expect_error(.check_autocor_family(f, stats::binomial(), "predict"),
               "prediction cannot proceed")
})

test_that(".check_autocor_family() passes everything brms can predict", {
  ar_plain <- brms::bf(y ~ x + ar(time = t))
  ar_cov   <- brms::bf(y ~ x + ar(time = t, cov = TRUE))
  no_ar    <- brms::bf(y ~ x)

  # Gaussian and student have natural residuals: cov = FALSE is fine. This is
  # the path every case study in the paper uses, so over-firing here would be
  # worse than the bug.
  expect_silent(.check_autocor_family(ar_plain, stats::gaussian(), "fit"))
  expect_silent(.check_autocor_family(ar_plain, brms::student(), "fit"))
  expect_silent(.check_autocor_family(ar_plain, NULL, "fit"))
  # Non-Gaussian is fine once cov = TRUE.
  expect_silent(.check_autocor_family(ar_cov, stats::poisson(), "fit"))
  # No autocorrelation term: the family is irrelevant.
  expect_silent(.check_autocor_family(no_ar, stats::poisson(), "fit"))
  expect_silent(.check_autocor_family(no_ar, stats::binomial(), "predict"))
})

test_that(".check_autocor_family() warns rather than errors when cov is unknown", {
  # A cov argument the constructors cannot evaluate must not hard-block a fit;
  # NA means "could not tell", and blocking on it would be worse than warning.
  fake <- structure(list(), class = "formula")
  expect_true(is.na(.autocor_cov(fake)))
})

# ── et_fit(): the gate fires before Stan, not after ────────────────────────
test_that("et_fit() refuses non-Gaussian ar(cov = FALSE) without compiling", {
  set.seed(1)
  d <- data.frame(y = rpois(30, 4), x = rnorm(30), t = 1:30)
  t0 <- Sys.time()
  expect_error(
    et_fit(y ~ x + ar(time = t, p = 1), data = d, family = stats::poisson(),
           chains = 1L, iter = 200L, silent = 2L, refresh = 0),
    "requires cov = TRUE"
  )
  # The whole point is failing fast. Stan compilation is tens of seconds at
  # best, so anything under a couple of seconds proves nothing was built.
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 5)
})

test_that("et_fit() accepts the family as a function or a string", {
  set.seed(1)
  d <- data.frame(y = rpois(30, 4), x = rnorm(30), t = 1:30)
  # family = poisson (unevaluated) and family = "poisson" are both legal brms
  # input, so the gate must normalise them rather than silently passing them.
  expect_error(
    et_fit(y ~ x + ar(time = t, p = 1), data = d, family = stats::poisson,
           chains = 1L, iter = 200L, silent = 2L, refresh = 0),
    "requires cov = TRUE"
  )
  expect_error(
    et_fit(y ~ x + ar(time = t, p = 1), data = d, family = "poisson",
           chains = 1L, iter = 200L, silent = 2L, refresh = 0),
    "requires cov = TRUE"
  )
})

# ── The cov = TRUE path: works, but says what it is handing back ───────────
skip_if_not_installed("brms")
if (requireNamespace("covr", quietly = TRUE) && covr::in_covr()) {
  Sys.setenv(NOT_CRAN = "true")
}
skip_on_cran()

test_that("non-Gaussian ar(cov = TRUE) warns that draws may leave the support", {
  set.seed(42)
  n  <- 50
  x  <- rnorm(n)
  e  <- as.numeric(stats::arima.sim(list(ar = 0.7), n, sd = 0.5))
  d  <- data.frame(y = rpois(n, exp(0.7 + 0.4 * x + e)), x = x, t = seq_len(n))
  tr <- d[1:38, ]
  vl <- d[39:n, ]

  fit <- et_fit(y ~ x + ar(time = t, p = 1, cov = TRUE), data = tr,
                family = stats::poisson(), chains = 2L, iter = 1000L,
                warmup = 500L, cores = 1L, seed = 1L, silent = 2L, refresh = 0)

  # With driver uncertainty the location shift is used and must announce the
  # support violation.
  expect_message(
    pred <- et_predict(fit, newdata = vl, env_noise = list(x = 0.5), seed = 1L),
    "NOT guaranteed to respect the response support"
  )

  # And the warning is not hypothetical: the draws really do leave the support
  # for a small-count Poisson.
  pd <- pred$predictive_draws
  expect_true(any(abs(pd - round(pd)) > 1e-9))

  # Without driver uncertainty there is no shift, so no warning and the draws
  # are genuine counts from brms.
  expect_silent(
    pred0 <- et_predict(fit, newdata = vl, env_noise = NULL, seed = 1L)
  )
  pd0 <- pred0$predictive_draws
  expect_true(all(abs(pd0 - round(pd0)) < 1e-9))
  expect_gte(min(pd0), 0)
})
