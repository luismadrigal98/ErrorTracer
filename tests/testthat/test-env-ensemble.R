# tests/testthat/test-env-ensemble.R
#
# Unit tests for the pathwise driver-ensemble path (T2): validation and the
# key property that a fanning-out ensemble makes the linear-predictor variance
# grow with lead time. No Stan fit needed.

# ── .validate_env_ensemble ───────────────────────────────────────────────────

test_that(".validate_env_ensemble accepts a well-formed named list", {
  nd  <- data.frame(x = 1:5, z = 6:10)
  ens <- list(x = matrix(rnorm(3 * 5), nrow = 3))   # 3 scenarios x 5 horizon
  out <- ErrorTracer:::.validate_env_ensemble(ens, c("x", "z"), nd)
  expect_type(out, "list")
  expect_equal(dim(out$x), c(3L, 5L))
})

test_that(".validate_env_ensemble errors on horizon mismatch", {
  nd  <- data.frame(x = 1:5)
  ens <- list(x = matrix(rnorm(3 * 4), nrow = 3))   # 4 != nrow(newdata)=5
  expect_error(ErrorTracer:::.validate_env_ensemble(ens, "x", nd),
               "ncol == nrow\\(newdata\\)")
})

test_that(".validate_env_ensemble errors on unnamed list and <2 scenarios", {
  nd <- data.frame(x = 1:5)
  expect_error(
    ErrorTracer:::.validate_env_ensemble(list(matrix(0, 3, 5)), "x", nd),
    "NAMED list")
  expect_error(
    ErrorTracer:::.validate_env_ensemble(list(x = matrix(0, 1, 5)), "x", nd),
    "at least 2 scenarios")
})

test_that(".validate_env_ensemble drops unknown predictors with a warning", {
  nd  <- data.frame(x = 1:5)
  ens <- list(x = matrix(0, 3, 5), WIND = matrix(0, 3, 5))
  # .et_warn() emits via message(), not warning()
  expect_message(
    out <- ErrorTracer:::.validate_env_ensemble(ens, "x", nd),
    "not in model")
  out <- suppressMessages(ErrorTracer:::.validate_env_ensemble(ens, "x", nd))
  expect_equal(names(out), "x")
})

test_that(".validate_env_ensemble errors when scenario counts differ", {
  nd  <- data.frame(x = 1:5, z = 6:10)
  ens <- list(x = matrix(0, 3, 5), z = matrix(0, 4, 5))
  expect_error(ErrorTracer:::.validate_env_ensemble(ens, c("x", "z"), nd),
               "same number of scenarios")
})

# ── .compute_lp_ensemble: variance grows with lead time ──────────────────────

test_that(".compute_lp_ensemble variance accumulates with a fanning ensemble", {
  set.seed(3)
  H <- 12L; M <- 50L; n_perturb <- 400L
  draws_mat <- cbind(b_Intercept = rnorm(n_perturb, 0, 0.05),
                     b_x          = rnorm(n_perturb, 1, 0.05))  # beta_x ~ 1
  x_mean <- seq(0, 1, length.out = H)
  # Random-walk scenarios: spread grows with lead time
  ens_x  <- t(replicate(M, x_mean + cumsum(rnorm(H, 0, 0.2))))
  nd     <- data.frame(x = x_mean)
  sid    <- rep_len(seq_len(M), n_perturb)

  lp <- ErrorTracer:::.compute_lp_ensemble(
    draws_mat = draws_mat, newdata = nd, pred_names = "x",
    env_ensemble = list(x = ens_x), scenario_ids = sid)

  expect_equal(dim(lp), c(n_perturb, H))
  v <- apply(lp, 2, var)
  # Variance at the far horizon should greatly exceed the near horizon
  expect_gt(v[H], 3 * v[1])
  # And it should be (roughly) monotonically increasing
  expect_gt(cor(seq_len(H), v), 0.9)
})

test_that(".compute_lp_ensemble holds non-ensemble predictors fixed", {
  set.seed(4)
  H <- 5L; M <- 10L; n_perturb <- 20L
  draws_mat <- cbind(b_Intercept = rep(0, n_perturb),
                     b_x = rep(1, n_perturb), b_z = rep(1, n_perturb))
  nd  <- data.frame(x = rep(0, H), z = 1:H)        # z fixed, x from ensemble
  ens <- list(x = matrix(2, M, H))                 # every scenario x = 2
  sid <- rep_len(seq_len(M), n_perturb)
  lp  <- ErrorTracer:::.compute_lp_ensemble(
    draws_mat = draws_mat, newdata = nd, pred_names = c("x", "z"),
    env_ensemble = ens, scenario_ids = sid)
  # lp = 0 + 2 (x from ensemble) + z (fixed) = 2 + (1:H)
  expect_equal(as.numeric(lp[1, ]), 2 + (1:H))
})
