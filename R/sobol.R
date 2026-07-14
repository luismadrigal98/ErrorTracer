# R/sobol.R -- et_sobol(): optional global (variance-based) sensitivity.
#
# The default decomposition (decompose_uncertainty) is a fast, interpretable
# ADDITIVE budget. It is exact under a linear link with independent factors,
# but -- as a reviewer noted -- a sequential additive split can be sensitive to
# interactions between the parameter and driver factors (e.g. under a non-linear
# inverse link, or with correlated drivers). et_sobol() provides the general,
# order-independent alternative: a variance-based (Sobol) decomposition of the
# predictive MEAN into first-order parameter and environmental contributions
# PLUS their interaction, using the Saltelli (2010) / Jansen estimators. The
# residual (within-draw family variance) is reported alongside so the full
# predictive variance is Var(mu) + residual, with Var(mu) split by Sobol.

#' Global (Sobol) variance-based sensitivity of the predictive mean
#'
#' An optional, order-independent alternative to the additive budget from
#' \code{\link{decompose_uncertainty}}.  Decomposes the variance of the
#' response-scale predictive mean \eqn{\mu = g^{-1}(\eta)} into first-order
#' \strong{parameter} and \strong{environmental} (driver) Sobol indices plus
#' their \strong{interaction}, by independently resampling the parameter draws
#' and the predictor perturbations (Saltelli 2010 first-order; Jansen
#' total-order).  Unlike the sequential additive split, this is insensitive to
#' the order in which factors are added and makes any parameter\eqn{\times}driver
#' interaction explicit --- the general case a non-linear link or correlated
#' drivers can produce.
#'
#' @param predictions An \code{et_prediction} that was produced with non-zero
#'   \code{env_noise} (the environmental factor needs a perturbation scale).
#' @param n_sobol Integer.  Base sample size for the Sobol estimator (default
#'   1024).  Total model evaluations per observation are \code{4 * n_sobol}.
#' @param seed Optional integer for reproducible sampling.
#' @param ... Unused.
#' @return A \code{data.frame} with one row per forecast observation:
#'   \describe{
#'     \item{obs_id}{Observation index.}
#'     \item{S_param, S_env}{First-order Sobol indices (fraction of the
#'       predictive-mean variance from parameters / drivers alone).}
#'     \item{ST_param, ST_env}{Total-order indices (including interactions).}
#'     \item{interaction}{\code{1 - S_param - S_env}: the share of
#'       predictive-mean variance attributable to the parameter\eqn{\times}driver
#'       interaction.}
#'     \item{var_mu}{Variance of the predictive mean (the quantity being
#'       decomposed).}
#'     \item{residual_var}{Mean within-draw family variance (the irreducible
#'       part; \code{var_mu + residual_var} approximates the total predictive
#'       variance).}
#'   }
#' @seealso \code{\link{decompose_uncertainty}} (the default additive budget).
#' @export
et_sobol <- function(predictions, n_sobol = 1024L, seed = NULL, ...) {
  UseMethod("et_sobol")
}

#' @export
et_sobol.et_prediction <- function(predictions, n_sobol = 1024L, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  n_sobol <- max(64L, as.integer(n_sobol))

  model   <- predictions$model
  fit     <- model$fit
  newdata <- predictions$newdata
  family  <- fit$family

  pred_names <- .brms_pred_names(fit)
  noise_sds  <- .resolve_env_noise(predictions$env_noise, pred_names, newdata)
  if (all(vapply(noise_sds, function(v) all(v == 0, na.rm = TRUE), logical(1)))) {
    stop("et_sobol() needs a non-zero env_noise: the environmental factor has ",
         "no perturbation scale to vary. Re-run et_predict() with env_noise.")
  }

  draws_mat   <- .brms_draws_matrix(fit)
  S           <- nrow(draws_mat)
  beta_cols   <- paste0("b_", pred_names)
  avail_betas <- intersect(beta_cols, colnames(draws_mat))
  avail_preds <- sub("^b_", "", avail_betas)
  has_int     <- "b_Intercept" %in% colnames(draws_mat)
  disp_draws  <- .extract_disp_draws(draws_mat, S)

  X   <- as.matrix(newdata[, avail_preds, drop = FALSE])   # n x p
  n_obs <- nrow(X); p <- length(avail_preds)
  N   <- n_sobol

  # Two independent parameter samples (A, B) shared across observations.
  iA <- sample.int(S, N, replace = TRUE)
  iB <- sample.int(S, N, replace = TRUE)
  betasA <- draws_mat[iA, avail_betas, drop = FALSE]
  betasB <- draws_mat[iB, avail_betas, drop = FALSE]
  intA <- if (has_int) draws_mat[iA, "b_Intercept"] else rep(0, N)
  intB <- if (has_int) draws_mat[iB, "b_Intercept"] else rep(0, N)

  # Response-scale mean for a parameter set + additive predictor perturbation.
  mu_of <- function(betas, intc, eps, xi) {
    xmat <- sweep(eps, 2, xi, "+")            # (x_i + eps): N x p
    lp   <- intc + rowSums(xmat * betas)
    as.numeric(.apply_inv_link(matrix(lp, ncol = 1L), family))
  }

  rows <- vector("list", n_obs)
  for (i in seq_len(n_obs)) {
    xi   <- X[i, ]
    sd_i <- vapply(avail_preds, function(pp) {
      v <- noise_sds[[pp]]; if (is.null(v)) 0 else v[min(i, length(v))]
    }, numeric(1))
    eA <- matrix(stats::rnorm(N * p), N, p) * matrix(sd_i, N, p, byrow = TRUE)
    eB <- matrix(stats::rnorm(N * p), N, p) * matrix(sd_i, N, p, byrow = TRUE)

    fA  <- mu_of(betasA, intA, eA, xi)        # (theta_A, eps_A)
    fB  <- mu_of(betasB, intB, eB, xi)        # (theta_B, eps_B)
    fCp <- mu_of(betasB, intB, eA, xi)        # param from B, env from A
    fCe <- mu_of(betasA, intA, eB, xi)        # param from A, env from B

    # Centre the outputs by the grand mean before the estimators. The Saltelli
    # first-order estimator mean(fB * (f_ABi - fA)) is inflated by a large
    # output mean (e.g. a large E[b]); centring is mean-invariant in theory and
    # sharply reduces the Monte Carlo variance in practice.
    f0  <- mean(c(fA, fB))
    fA  <- fA - f0; fB <- fB - f0; fCp <- fCp - f0; fCe <- fCe - f0

    V <- stats::var(c(fA, fB))
    if (!is.finite(V) || V < 1e-12) {
      S_param <- S_env <- ST_param <- ST_env <- interaction <- NA_real_
    } else {
      S_param <- mean(fB * (fCp - fA)) / V             # Saltelli 2010
      S_env   <- mean(fB * (fCe - fA)) / V
      ST_param <- 0.5 * mean((fA - fCp)^2) / V         # Jansen
      ST_env   <- 0.5 * mean((fA - fCe)^2) / V
      interaction <- 1 - S_param - S_env
    }

    rows[[i]] <- data.frame(
      obs_id = i,
      S_param = S_param, S_env = S_env,
      ST_param = ST_param, ST_env = ST_env,
      interaction = interaction,
      var_mu = V,
      residual_var = NA_real_,   # filled below (family variance)
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)

  # Residual (within-draw family variance), consistent with the additive budget.
  mu_draws <- if (.is_gaussian_identity(family)) {
    predictions$posterior_linpred
  } else {
    .apply_inv_link(predictions$posterior_linpred, family)
  }
  out$residual_var <- .family_residual_var(family, mu_draws, disp_draws)

  structure(out, class = c("et_sobol", "data.frame"), n_sobol = N)
}

#' @export
et_sobol.default <- function(predictions, ...) {
  stop("et_sobol() expects an et_prediction object (with non-zero env_noise).")
}

#' @export
print.et_sobol <- function(x, ...) {
  cat("ErrorTracer global (Sobol) sensitivity of the predictive mean\n")
  cat("  Observations :", nrow(x), "\n")
  cat("  Base samples :", attr(x, "n_sobol"), "\n")
  cat(sprintf("  Mean first-order  : param = %.3f   env = %.3f\n",
              mean(x$S_param, na.rm = TRUE), mean(x$S_env, na.rm = TRUE)))
  cat(sprintf("  Mean total-order  : param = %.3f   env = %.3f\n",
              mean(x$ST_param, na.rm = TRUE), mean(x$ST_env, na.rm = TRUE)))
  cat(sprintf("  Mean interaction  : %.3f  (param x env share of mean variance)\n",
              mean(x$interaction, na.rm = TRUE)))
  invisible(x)
}
