# tests/testthat/test-decompose.R

# Build a minimal mock et_prediction object to test decompose_uncertainty
# without requiring a real Stan/brms fit.
.mock_et_prediction <- function(n_obs = 5, n_draws = 100) {
  set.seed(99)
  pp <- matrix(rnorm(n_draws * n_obs, mean = 0.5), nrow = n_draws)
  lp <- matrix(rnorm(n_draws * n_obs, mean = 0.5, sd = 0.3), nrow = n_draws)
  # Perturbed LP — slightly higher variance than lp
  lp_p <- lp + matrix(rnorm(n_draws * n_obs, sd = 0.1), nrow = n_draws)
  sigma_draws <- abs(rnorm(n_draws, mean = 0.5, sd = 0.1))

  decomp <- ErrorTracer:::.decompose_from_arrays(pp, lp, lp_p, sigma_draws)

  ci_df <- ErrorTracer:::.compute_ci(pp, c(0.5, 0.9, 0.95))

  structure(
    list(
      posterior_predict  = pp,
      posterior_linpred  = lp,
      lp_perturbed       = lp_p,
      sigma_draws        = sigma_draws,
      credible_intervals = ci_df,
      decomposition      = decomp,
      newdata            = data.frame(obs_id = seq_len(n_obs)),
      model              = NULL,
      env_noise          = NULL,
      n_draws            = n_draws
    ),
    class = "et_prediction"
  )
}

test_that("decompose_uncertainty returns correct columns", {
  pred   <- .mock_et_prediction()
  decomp <- decompose_uncertainty(pred)

  expect_s3_class(decomp, "data.frame")
  expect_true(all(c("obs_id", "param_var", "env_var",
                    "residual_var", "total_var") %in% colnames(decomp)))
})

test_that("decompose_uncertainty has non-negative variance components", {
  pred   <- .mock_et_prediction(n_obs = 10)
  decomp <- decompose_uncertainty(pred)

  expect_true(all(decomp$param_var    >= 0))
  expect_true(all(decomp$env_var      >= 0))
  expect_true(all(decomp$residual_var >= 0))
  expect_true(all(decomp$total_var    >= 0))
})

test_that("decompose_uncertainty row count matches n_obs", {
  n_obs <- 7
  pred  <- .mock_et_prediction(n_obs = n_obs)
  decomp <- decompose_uncertainty(pred)
  expect_equal(nrow(decomp), n_obs)
})

test_that("decompose_uncertainty.default raises an error", {
  expect_error(decompose_uncertainty(list(a = 1)), "et_prediction")
})

test_that("param_var < total_var on average (residual inflates total)", {
  pred   <- .mock_et_prediction(n_obs = 20, n_draws = 500)
  decomp <- decompose_uncertainty(pred)
  # total_var includes residual noise, so should generally exceed param_var
  expect_true(mean(decomp$total_var) > mean(decomp$param_var))
})

test_that("residual_var is constant across observations (scalar sigma)", {
  pred   <- .mock_et_prediction(n_obs = 6)
  decomp <- decompose_uncertainty(pred)
  # residual_var is median(sigma^2) — same for all obs
  expect_equal(length(unique(round(decomp$residual_var, 10))), 1L)
})
