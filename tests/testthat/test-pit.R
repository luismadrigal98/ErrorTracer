# tests/testthat/test-pit.R
#
# Unit tests for the PIT calibration diagnostics (T4): et_pit() core logic,
# the et_prediction method, and et_plot_pit(). No Stan fit needed.

test_that(".compute_pit computes the empirical predictive CDF", {
  # Each obs's predictive column is 1..100; non-randomized PIT = P(Y <= y).
  draws <- matrix(rep(1:100, 3), ncol = 3)
  y     <- c(25, 50, 75)
  pit   <- ErrorTracer:::.compute_pit(draws, y, randomize = FALSE)
  expect_equal(pit, c(0.25, 0.50, 0.75), tolerance = 1e-8)
})

test_that(".compute_pit returns NA for NA observations", {
  draws <- matrix(rnorm(100 * 2), ncol = 2)
  pit   <- ErrorTracer:::.compute_pit(draws, c(NA, 0), randomize = FALSE)
  expect_true(is.na(pit[1]))
  expect_false(is.na(pit[2]))
})

test_that("et_pit returns pit in [0,1]; et_plot_pit returns a ggplot", {
  set.seed(1)
  S <- 400L; n <- 20L
  pp   <- matrix(rnorm(S * n, mean = 0.5), nrow = S)
  pred <- structure(list(
    predictive_draws  = pp,
    posterior_predict = pp,
    newdata           = data.frame(x = rnorm(n)),
    model             = NULL
  ), class = "et_prediction")
  obs <- data.frame(y = rnorm(n, 0.5))

  pit <- et_pit(pred, observed = obs, response_col = "y")
  expect_true(all(pit$pit >= 0 & pit$pit <= 1))
  expect_equal(nrow(pit), n)
  expect_named(pit, c("obs_id", "pit"))

  expect_s3_class(et_plot_pit(pit), "ggplot")
  expect_s3_class(et_plot_pit(runif(50)), "ggplot")   # bare numeric vector
})

test_that("et_pit falls back to posterior_predict when predictive_draws absent", {
  set.seed(3)
  S <- 200L; n <- 8L
  pp   <- matrix(rnorm(S * n), nrow = S)
  pred <- structure(list(
    posterior_predict = pp,      # no predictive_draws (older object)
    newdata           = data.frame(x = rnorm(n)),
    model             = NULL
  ), class = "et_prediction")
  pit <- et_pit(pred, observed = data.frame(y = rnorm(n)), response_col = "y")
  expect_equal(nrow(pit), n)
  expect_true(all(is.finite(pit$pit)))
})

test_that("randomized PIT is ~uniform for a discrete (Poisson) predictive", {
  set.seed(2)
  S <- 2000L; n <- 300L; lambda <- 5
  pp   <- matrix(rpois(S * n, lambda), nrow = S)
  y    <- rpois(n, lambda)
  pred <- structure(list(predictive_draws = pp,
                         newdata = data.frame(z = seq_len(n)), model = NULL),
                    class = "et_prediction")

  pit_r <- et_pit(pred, observed = data.frame(y = y), response_col = "y",
                  randomize = TRUE)
  ks <- suppressWarnings(stats::ks.test(pit_r$pit, "punif"))
  expect_gt(ks$p.value, 0.01)          # randomized PIT ~ uniform
  expect_true(all(pit_r$pit >= 0 & pit_r$pit <= 1))
})
