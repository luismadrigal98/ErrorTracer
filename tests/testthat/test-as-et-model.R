# tests/testthat/test-as-et-model.R
#
# as_et_model() and the brmsfit methods added in 1.4.2: a model fitted outside
# the package enters the pipeline without being refitted.
#
# The abstract of the accompanying paper claims that any existing brms fit can
# be supplied directly. Until 1.4.2 it could not: et_predict() dispatched on
# et_model / et_model_list only, so a brmsfit raised "no applicable method".
#
# The error paths need no Stan and run everywhere; the round-trip tests fit a
# small model and are skipped on CRAN with the rest of the Stan-backed suite.

skip_if_not_installed("brms")
options(brms.backend = "rstan")

# ── Error paths: no fit required ─────────────────────────────────────────────

test_that("a bare stanfit is refused, with the brm(empty = TRUE) route", {
  sf <- structure(list(), class = "stanfit")
  expect_error(as_et_model(sf), "empty = TRUE", fixed = TRUE)
  expect_error(as_et_model(sf), "carries none of them")
  # The same message, whichever Stan object is handed over.
  expect_error(as_et_model(structure(list(), class = "CmdStanMCMC")),
               "empty = TRUE", fixed = TRUE)
  expect_error(as_et_model(structure(list(), class = "stanmodel")),
               "empty = TRUE", fixed = TRUE)
})

test_that("an object that is not a brms fit is refused by class", {
  d <- data.frame(y = rnorm(10), x = rnorm(10))
  expect_error(as_et_model(stats::lm(y ~ x, data = d)), "needs a brmsfit")
  expect_error(as_et_model(42), "needs a brmsfit")
})

test_that("a multivariate fit is refused before anything else is read", {
  fake <- structure(
    list(formula = structure(list(), class = c("mvbrmsformula", "bform"))),
    class = "brmsfit")
  expect_error(as_et_model(fake), "[Mm]ultivariate")
})

test_that("a fit carrying no training data is refused", {
  fake <- structure(list(formula = y ~ x, data = NULL), class = "brmsfit")
  expect_error(as_et_model(fake), "no training data")
  fake$data <- data.frame()
  expect_error(as_et_model(fake), "no training data")
})

# ── Round trip through a real fit ────────────────────────────────────────────

.fit_env <- new.env(parent = emptyenv())

.external_fit <- function() {
  if (is.null(.fit_env$fit)) {
    set.seed(11L)
    x <- rnorm(40)
    d <- data.frame(y = rnorm(40, 1 + 2 * x, 0.5), x = x)
    .fit_env$data <- d
    .fit_env$fit  <- brms::brm(
      y ~ x, data = d, chains = 1L, iter = 600L, warmup = 300L,
      seed = 4242L, refresh = 0L, silent = 2L
    )
  }
  .fit_env$fit
}

test_that("as_et_model() recovers formula, data and sampling seed", {
  skip_on_cran()
  m <- as_et_model(.external_fit())

  expect_s3_class(m, "et_model")
  expect_setequal(names(m),
                  c("fit", "prior_spec", "formula", "data", "config"))
  expect_s3_class(m$fit, "brmsfit")
  expect_null(m$prior_spec)
  expect_equal(nrow(m$data), 40L)
  expect_equal(m$config$seed, 4242L)
  expect_equal(m$config$chains, 1L)
})

test_that("the stored formula is the plain one, so the response is findable", {
  skip_on_cran()
  m <- as_et_model(.external_fit())
  # brms stores a brmsformula, and stats::update() on one yields nothing
  # usable -- which is exactly how et_calibrate() locates the response column.
  expect_false(inherits(m$formula, "brmsformula"))
  expect_equal(all.vars(stats::update(m$formula, . ~ 0)), "y")
  expect_length(all.vars(stats::update(.external_fit()$formula, . ~ 0)), 0L)
})

test_that("et_predict() accepts a brmsfit and matches the wrapped model", {
  skip_on_cran()
  fit <- .external_fit()
  nd  <- head(.fit_env$data, 6L)

  direct  <- et_predict(fit, newdata = nd, env_noise = list(x = 0.1),
                        n_draws = 200L)
  wrapped <- et_predict(as_et_model(fit), newdata = nd,
                        env_noise = list(x = 0.1), n_draws = 200L)

  expect_s3_class(direct, "et_prediction")
  # The seed is recovered from the fit, so both paths draw the same stream.
  expect_equal(direct$predictive_draws, wrapped$predictive_draws)
  expect_equal(direct$decomposition$total_var,
               wrapped$decomposition$total_var)
})

test_that("the pipeline runs end to end from an external fit", {
  skip_on_cran()
  fit <- .external_fit()
  nd  <- head(.fit_env$data, 8L)

  pred <- et_predict(fit, newdata = nd, env_noise = list(x = 0.1),
                     n_draws = 300L)
  dec  <- decompose_uncertainty(pred)
  expect_true(all(c("param_var", "env_var", "residual_var", "total_var")
                  %in% colnames(dec)))
  expect_true(all(dec$total_var > 0))
  expect_budget_closes(dec, pred)

  sl <- shelf_life(pred, response_scale = range(.fit_env$data$y))
  expect_s3_class(sl, "et_shelf_life")
  expect_true(all(sl$ratio > 0))

  cal <- et_calibrate(pred, observed = nd)
  expect_true(all(cal$observed_coverage >= 0 & cal$observed_coverage <= 1))
})

test_that("seed is overridable and a bogus prior_spec is dropped", {
  skip_on_cran()
  fit <- .external_fit()
  expect_equal(as_et_model(fit, seed = 7L)$config$seed, 7L)
  expect_null(as_et_model(fit, prior_spec = list(method = "nonsense"))$prior_spec)
  # An explicit NA seed is how et_predict() is told to use the ambient stream.
  expect_true(is.na(as_et_model(fit, seed = NA_integer_)$config$seed))
})

test_that("supplied data overrides the fit's model frame", {
  skip_on_cran()
  fit <- .external_fit()
  d2  <- cbind(.fit_env$data, extra = 1)
  expect_true("extra" %in% colnames(as_et_model(fit, data = d2)$data))
})

test_that("et_diagnose() accepts a brmsfit", {
  skip_on_cran()
  fit <- .external_fit()
  dg  <- et_diagnose(fit, loo = FALSE)
  expect_equal(dg, et_diagnose(as_et_model(fit), loo = FALSE))
  expect_true(!is.null(dg$convergence))
})

test_that("print() on a wrapped model reports the recovered settings", {
  skip_on_cran()
  out <- capture.output(print(as_et_model(.external_fit())))
  expect_true(any(grepl("et_model", out)))
  expect_true(any(grepl("brms defaults", out)))
  expect_false(any(grepl("Chains  :\\s*Iter", out)))
})
