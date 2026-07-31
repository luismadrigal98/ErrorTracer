# Data-splitting prior extraction (ErrorTracer 1.3.1).
#
# R1 (MEE-26-05-477, major flaw 1): "Step 1 of the workflow (extracting
# coefficients, using them as priors) violates the fundamental separation
# between prior information and information entering through the likelihood."
#
# shrinkage = "zero" removed the reuse of the point estimate as a prior MEAN.
# What remained was the prior scale, and -- more consequentially -- variable
# SELECTION performed on the same rows that form the likelihood. These tests pin
# that et_priors_split() removes both by construction.

test_that("the two halves are disjoint and exhaustive", {
  set.seed(1)
  d <- data.frame(x1 = rnorm(60), x2 = rnorm(60))
  d$y <- 1 + 2 * d$x1 + rnorm(60)

  sp <- et_priors_split(d, function(df) lm(y ~ x1 + x2, data = df), seed = 42)

  expect_equal(sp$n_prior + sp$n_fit, nrow(d))
  # The substantive property: no observation informs both the prior and the
  # likelihood.
  expect_equal(nrow(merge(sp$prior_data, sp$fit_data)), 0)
})

test_that("the split is reproducible under a seed", {
  set.seed(1)
  d <- data.frame(x1 = rnorm(50)); d$y <- d$x1 + rnorm(50)
  f <- function(df) lm(y ~ x1, data = df)
  a <- et_priors_split(d, f, seed = 7)
  b <- et_priors_split(d, f, seed = 7)
  expect_identical(rownames(a$fit_data), rownames(b$fit_data))
})

test_that("a supplied seed does not leak into the caller's RNG stream", {
  set.seed(1)
  d <- data.frame(x1 = rnorm(50)); d$y <- d$x1 + rnorm(50)

  set.seed(99); before <- rnorm(3)
  set.seed(99)
  invisible(et_priors_split(d, function(df) lm(y ~ x1, data = df), seed = 7))
  after <- rnorm(3)

  # Without state restoration the internal set.seed(7) would reseed the session
  # and these would differ -- a silent reproducibility break for the caller.
  expect_equal(before, after)
})

test_that("prop controls the split size", {
  set.seed(1)
  d <- data.frame(x1 = rnorm(100)); d$y <- d$x1 + rnorm(100)
  sp <- et_priors_split(d, function(df) lm(y ~ x1, data = df),
                        prop = 0.3, seed = 1)
  expect_equal(sp$n_prior, 30)
  expect_equal(sp$n_fit, 70)
})

test_that("stratified splitting represents every level on both sides", {
  set.seed(1)
  d <- data.frame(g = rep(c("a", "b", "c"), each = 20), x1 = rnorm(60))
  d$y <- d$x1 + rnorm(60)
  sp <- et_priors_split(d, function(df) lm(y ~ x1, data = df),
                        strata = "g", seed = 3)
  expect_setequal(unique(sp$prior_data$g), c("a", "b", "c"))
  expect_setequal(unique(sp$fit_data$g),   c("a", "b", "c"))
})

test_that("priors are extracted from the prior half only", {
  set.seed(1)
  d <- data.frame(x1 = rnorm(60)); d$y <- 5 * d$x1 + rnorm(60)
  seen <- NULL
  f <- function(df) { seen <<- nrow(df); lm(y ~ x1, data = df) }
  sp <- et_priors_split(d, f, prop = 0.4, seed = 1)
  expect_equal(seen, 24)          # prior_fun never sees the likelihood rows
  expect_s3_class(sp$priors, "et_prior_spec")
})

test_that("degenerate splits error rather than silently returning nothing", {
  d <- data.frame(x1 = rnorm(10)); d$y <- d$x1 + rnorm(10)
  f <- function(df) lm(y ~ x1, data = df)
  expect_error(et_priors_split(d, f, prop = 0.99), "row\\(s\\) for the likelihood")
  expect_error(et_priors_split(d, f, prop = 1.5), "strictly between 0 and 1")
  expect_error(et_priors_split(d[1:3, ], f), "at least 4 rows")
})
