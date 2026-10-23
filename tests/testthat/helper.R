# tests/testthat/helper.R
# Shared setup for testthat tests — loaded automatically by testthat

# Suppress brms/Stan messages during tests
options(brms.backend = "rstan")

# ── Budget closure ──────────────────────────────────────────────────────────
# Since 1.4.0 total_var is the MEASURED variance of the env-inclusive posterior
# predictive draws, not the sum of the channels, so the channels close on it as
# a RESULT rather than by definition. The closure gap is therefore Monte Carlo
# error, and these helpers bound it by its own scale instead of asserting an
# identity to 1e-8 (which the old definitional total satisfied tautologically).
#
# Scale: a variance estimated from n draws has relative SE sqrt(2 / (n - 1)),
# and the gap combines three such estimates (Var(pp_env), param_var, and the
# env subtraction over n_perturb draws). Bound the MEDIAN relative gap at
# 3 * sqrt(2 / n_eff) so one noisy observation cannot fail the suite.
#
# NOTE: for non-Gaussian families there is additionally a small STRUCTURAL gap
# (measured at 1-2% of total_var on a Poisson/log fit), because the
# env-inclusive draws are built as a location shift
# mu_env + (y_pp - mu) rather than as a draw from the family at mu_env. That is
# well inside the bound below; see .inflate_env_predictive().
expect_budget_closes <- function(dec, pred = NULL, n_eff = NULL,
                                  factor = 3) {
  comp <- dec$param_var + dec$env_var + dec$residual_var
  if ("group_var" %in% colnames(dec))    comp <- comp + dec$group_var
  if ("temporal_var" %in% colnames(dec)) comp <- comp + dec$temporal_var

  # With a temporal channel the closure is exact by construction (temporal_var
  # IS the remainder), so assert the identity there.
  if ("temporal_var" %in% colnames(dec)) {
    testthat::expect_equal(dec$total_var, comp, tolerance = 1e-8)
    return(invisible(NULL))
  }

  if (is.null(n_eff)) {
    n_eff <- if (!is.null(pred)) {
      min(pred$n_draws, pred$n_perturb %||% pred$n_draws)
    } else {
      400L
    }
  }
  tol <- factor * sqrt(2 / max(n_eff - 1, 1))
  rel <- abs(dec$total_var - comp) / pmax(abs(dec$total_var),
                                          .Machine$double.eps)
  testthat::expect_lt(stats::median(rel), tol)
  # closure_gap must be exactly the shortfall it claims to report.
  if ("closure_gap" %in% colnames(dec)) {
    testthat::expect_equal(dec$closure_gap, dec$total_var - comp)
  }
  invisible(NULL)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
