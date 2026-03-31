# tests/testthat/test-shelf_life.R

# Build a minimal mock et_prediction for shelf life tests
.mock_pred_for_sl <- function(n_obs = 6, ci_widths = NULL) {
  set.seed(42)
  n_draws <- 200
  # Widths grow with obs_id (simulating increasing uncertainty over forecast horizon)
  if (is.null(ci_widths)) ci_widths <- seq(0.2, 1.4, length.out = n_obs)

  # Construct a pp matrix whose quantile widths approximate ci_widths
  pp <- matrix(NA_real_, n_draws, n_obs)
  for (j in seq_len(n_obs)) {
    pp[, j] <- stats::rnorm(n_draws, sd = ci_widths[j] / (2 * 1.645))
  }

  lp  <- pp + matrix(rnorm(n_draws * n_obs, sd = 0.01), n_draws, n_obs)
  lp_p <- lp
  sigma_draws <- rep(0.1, n_draws)

  decomp <- ErrorTracer:::.decompose_from_arrays(pp, lp, lp_p, sigma_draws)
  ci_df  <- ErrorTracer:::.compute_ci(pp, c(0.90, 0.95))

  nd <- data.frame(year = 2020 + seq_len(n_obs) - 1L)

  structure(
    list(
      posterior_predict  = pp,
      posterior_linpred  = lp,
      lp_perturbed       = lp_p,
      sigma_draws        = sigma_draws,
      credible_intervals = ci_df,
      decomposition      = decomp,
      newdata            = nd,
      model              = NULL
    ),
    class = "et_prediction"
  )
}

test_that("shelf_life returns et_shelf_life with correct columns", {
  pred <- .mock_pred_for_sl()
  sl   <- shelf_life(pred, plausible_range = c(-1, 1), ci_level = 0.90)

  expect_s3_class(sl, c("et_shelf_life", "data.frame"))
  expect_true(all(c("obs_id", "time", "ci_width", "plausible_range",
                    "ratio", "informative") %in% colnames(sl)))
})

test_that("shelf_life row count matches n_obs", {
  n <- 8
  pred <- .mock_pred_for_sl(n_obs = n)
  sl   <- shelf_life(pred, plausible_range = c(-2, 2))
  expect_equal(nrow(sl), n)
})

test_that("shelf_life ratio = ci_width / plausible_range", {
  pred <- .mock_pred_for_sl()
  sl   <- shelf_life(pred, plausible_range = c(-1, 1))
  expect_equal(sl$ratio, sl$ci_width / 2, tolerance = 1e-10)
})

test_that("shelf_life informative flag respects threshold", {
  pred <- .mock_pred_for_sl()
  sl   <- shelf_life(pred, plausible_range = c(-1, 1), threshold = 0.5)
  # informative == TRUE when ratio < 0.5
  expect_equal(sl$informative, sl$ratio < 0.5)
})

test_that("shelf_life uses time_col from newdata", {
  pred <- .mock_pred_for_sl(n_obs = 4)
  sl   <- shelf_life(pred, plausible_range = c(-1, 1), time_col = "year")
  expect_equal(sl$time, 2020:2023)
})

test_that("shelf_life errors if ci_level not in predictions", {
  pred <- .mock_pred_for_sl()
  expect_error(
    shelf_life(pred, plausible_range = c(-1, 1), ci_level = 0.99),
    "ci_level"
  )
})

test_that("shelf_life errors if plausible_range has length != 2", {
  pred <- .mock_pred_for_sl()
  expect_error(shelf_life(pred, plausible_range = c(1, 2, 3)), "length 2")
})

test_that("shelf_life errors if plausible_range min == max", {
  pred <- .mock_pred_for_sl()
  expect_error(shelf_life(pred, plausible_range = c(1, 1)), "equal")
})

test_that("shelf_life print method runs without error", {
  pred <- .mock_pred_for_sl()
  sl   <- shelf_life(pred, plausible_range = c(-1, 1))
  expect_output(print(sl), "ErrorTracer shelf life")
  expect_output(print(sl), "Informative")
})

test_that("shelf_life plausible_range column is constant (scalar diff)", {
  pred <- .mock_pred_for_sl(n_obs = 5)
  sl   <- shelf_life(pred, plausible_range = c(-2, 2))
  expect_true(all(sl$plausible_range == 4))
})
