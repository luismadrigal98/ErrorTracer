# tests/testthat/test-decompose-autocor.R
#
# Regression tests for the variance decomposition under a residual
# autocorrelation term. These pin the two properties that the internal referee
# review of 2026-07-27 found violated:
#
#   M1  env_var must recover the analytic environmental variance beta^2 * delta^2
#       under ar() just as it does under an iid model. It currently does not:
#       v_perturbed comes from .compute_lp_perturbed() (a hand-rolled linear
#       predictor with NO ar term) while v_param_sub comes from
#       posterior_linpred() (which, under brms's default ar(cov = FALSE), DOES
#       carry the ar contribution to the mean). The difference is therefore
#       V_env - V_AR, which goes negative and is clamped by pmax(., 0).
#
#   M2  temporal_var must be > 0 for a series with strong residual
#       autocorrelation. It is currently 0 at every lead when phi ~ 0.8, because
#       param_var (computed from posterior_linpred on newdata alone, i.e. a
#       cold start) has already absorbed the ar-error variance, so
#       pp_var - (param_var + residual_var) <= 0.
#
# Note why the existing suite missed both: .mock_et_prediction() in
# test-decompose.R builds mu_perturbed as `lp + noise`, i.e. it *assumes* the
# two arrays are on the same footing. Only a real ar() fit exercises the
# asymmetry.
#
# Reproduction scripts: tests/manual/check1_ar_accumulation.R,
#                       tests/manual/check2_env_recovery.R
#
# Runtime: ~2-4 min on first call (Stan compilation); <30 s on reruns.

skip_if_not_installed("brms")
if (requireNamespace("covr", quietly = TRUE) && covr::in_covr()) {
  Sys.setenv(NOT_CRAN = "true")
}
skip_on_cran()

options(brms.backend = "rstan")

# ── Shared fixture: an AR(1) series with a known regression slope ────────────
# y_t = 1 + 2*x_t + e_t,  e_t = phi e_{t-1} + eta_t,  phi = 0.9, sd(eta) = 0.5
# Forecast x values are iid (NOT trending), so any growth in param_var with
# lead time cannot be extrapolation — it can only be ar error.

set.seed(42L)
.n <- 60L
.H <- 15L
.x <- rnorm(.n)
.phi <- 0.9
.sig <- 0.5
.e <- numeric(.n)
.e[1] <- rnorm(1, 0, .sig / sqrt(1 - .phi^2))
for (t in 2:.n) .e[t] <- .phi * .e[t - 1] + rnorm(1, 0, .sig)
.y <- 1 + 2 * .x + .e

.train <- data.frame(y = .y, x = .x, t = seq_len(.n))
.newd  <- data.frame(y = NA_real_, x = rnorm(.H), t = .n + seq_len(.H))

.delta <- 0.2   # env_noise SD on x

.fit_ar <- et_fit(y ~ x + ar(time = t, p = 1), data = .train,
                  chains = 2L, iter = 1500L, warmup = 750L, cores = 2L,
                  refresh = 0, seed = 1L, priors = NULL, silent = 2L)

.pred_ar <- et_predict(.fit_ar, newdata = .newd,
                       env_noise = list(x = .delta),
                       n_draws = 1500L, n_perturb = 500L)

.dec_ar   <- decompose_uncertainty(.pred_ar)
.beta_ar  <- mean(as.matrix(.fit_ar$fit, variable = "b_x")[, 1])
.v_env_true <- .beta_ar^2 * .delta^2     # analytic target, constant across leads


# ── M1: env_var must survive an ar() term ───────────────────────────────────

test_that("env_var recovers beta^2*delta^2 under an ar() term (M1)", {
  # The environmental variance does not depend on lead time here: env_noise is
  # constant and the forecast x values are exchangeable. Allow a generous
  # tolerance for Monte Carlo error and the Var(beta)*delta^2 term.
  expect_equal(mean(.dec_ar$env_var), .v_env_true, tolerance = 0.35)
})

test_that("env_var is not systematically floored to zero under ar() (M1)", {
  # The pmax(., 0) floor exists to absorb Monte Carlo jitter around a true zero.
  # With a genuine, sizeable env signal it should essentially never fire.
  n_floored <- sum(.dec_ar$env_var == 0)
  expect_lt(n_floored, 0.2 * nrow(.dec_ar))
})

test_that("env_var does not decay with lead time under ar() (M1)", {
  # The failure signature is env_var surviving at lead 1 and collapsing as the
  # ar error grows. Compare the first third of the horizon to the last third.
  k <- floor(nrow(.dec_ar) / 3)
  early <- mean(.dec_ar$env_var[seq_len(k)])
  late  <- mean(.dec_ar$env_var[(nrow(.dec_ar) - k + 1):nrow(.dec_ar)])
  expect_gt(late, 0.5 * early)
})


# ── M2: param_var / temporal_var must not be conflated ──────────────────────

test_that("temporal_var is positive for a strongly autocorrelated series (M2)", {
  skip_if(is.null(.dec_ar$temporal_var), "no temporal_var column returned")
  expect_gt(mean(.dec_ar$temporal_var), 0)
})

test_that("param_var does not absorb the ar error (M2)", {
  # Forecast x is iid, so regression parameter variance is flat in lead time.
  # Any strong monotone growth in param_var is ar error leaking into the
  # 'parameter' channel.
  fit <- lm(.dec_ar$param_var ~ seq_len(nrow(.dec_ar)))
  slope <- unname(coef(fit)[2])
  # Growth of more than ~2% of the mean param_var per lead step indicates leakage.
  expect_lt(slope, 0.02 * mean(.dec_ar$param_var))
})


# ── Control: the iid path must keep working ─────────────────────────────────

test_that("env_var recovers beta^2*delta^2 under an iid model (control)", {
  fit_iid <- et_fit(y ~ x, data = .train,
                    chains = 2L, iter = 1500L, warmup = 750L, cores = 2L,
                    refresh = 0, seed = 1L, priors = NULL, silent = 2L)
  pred_iid <- et_predict(fit_iid, newdata = .newd,
                         env_noise = list(x = .delta),
                         n_draws = 1500L, n_perturb = 500L)
  dec_iid  <- decompose_uncertainty(pred_iid)
  beta_iid <- mean(as.matrix(fit_iid$fit, variable = "b_x")[, 1])

  expect_equal(mean(dec_iid$env_var), beta_iid^2 * .delta^2, tolerance = 0.35)
  expect_equal(sum(dec_iid$env_var == 0), 0L)
})


# ── Budget consistency must hold regardless ─────────────────────────────────

test_that("components still reconcile to total_var exactly", {
  comp <- .dec_ar$param_var + .dec_ar$env_var + .dec_ar$residual_var
  if (!is.null(.dec_ar$temporal_var)) comp <- comp + .dec_ar$temporal_var
  expect_equal(comp, .dec_ar$total_var, tolerance = 1e-10)
})
