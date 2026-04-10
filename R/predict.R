# R/predict.R — et_predict(): posterior prediction with uncertainty propagation

# ============================================================
# S3 generic
# ============================================================

#' Posterior prediction with uncertainty decomposition
#'
#' Generates posterior predictive draws for new observations, propagates
#' environmental measurement uncertainty through the model, and computes
#' credible intervals.  The resulting \code{et_prediction} object is the
#' input to \code{\link{decompose_uncertainty}}, \code{\link{shelf_life}},
#' and the plotting functions.
#'
#' @param model An \code{et_model} or \code{et_model_list} object from
#'   \code{\link{et_fit}}.
#' @param newdata A \code{data.frame} containing the predictor columns
#'   named in the model formula.  For grouped models, must also contain
#'   the grouping column.
#' @param env_noise Environmental measurement / prediction uncertainty.
#'   Can be:
#'   \itemize{
#'     \item A named \code{list} or named numeric vector: per-predictor
#'       absolute noise SDs, e.g.\
#'       \code{list(Tmean = 0.5, PPT = 10)}.
#'     \item A single numeric: applied as a fraction of each predictor's
#'       empirical SD in \code{newdata} (e.g.\ \code{0.1} means 10\% noise).
#'     \item \code{NULL} (default): no environmental noise.
#'   }
#' @param n_draws Integer.  Number of posterior draws to use (default 2000;
#'   capped at the number of draws available in the fit).
#' @param ci_levels Numeric vector.  Credible interval levels to compute
#'   (default \code{c(0.5, 0.8, 0.9, 0.95)}).
#' @param n_perturb Integer.  Number of draws used for the environmental
#'   perturbation step (default \code{min(500, n_draws)}).  Reducing this
#'   speeds up computation.
#' @param ... Passed to methods.
#'
#' @return An \code{et_prediction} object (list) containing:
#' \describe{
#'   \item{\code{posterior_predict}}{Matrix \code{[n_draws x n_obs]}: full
#'     posterior predictive draws (parameter + residual uncertainty).}
#'   \item{\code{posterior_linpred}}{Matrix \code{[n_draws x n_obs]}: linear
#'     predictor draws (parameter uncertainty only).}
#'   \item{\code{lp_perturbed}}{Matrix \code{[n_perturb x n_obs]}: linear
#'     predictor computed on environmentally perturbed inputs.}
#'   \item{\code{sigma_draws}}{Numeric vector: posterior draws of sigma.}
#'   \item{\code{credible_intervals}}{data.frame with columns
#'     \code{row_id, ci_level, lower, median, upper, width}.}
#'   \item{\code{decomposition}}{data.frame from
#'     \code{decompose_uncertainty()}: \code{param_var, env_var,
#'     residual_var, total_var}.}
#'   \item{\code{newdata}}{The input \code{newdata}.}
#'   \item{\code{model}}{Reference to the \code{et_model} used.}
#' }
#'
#' @seealso \code{\link{decompose_uncertainty}}, \code{\link{shelf_life}},
#'   \code{\link{et_calibrate}}
#' @export
et_predict <- function(model, newdata, env_noise = NULL,
                        n_draws = 2000L, ci_levels = c(0.5, 0.8, 0.9, 0.95),
                        n_perturb = NULL, ...) {
  UseMethod("et_predict")
}

# ============================================================
# Method: et_model (single model)
# ============================================================

#' @export
et_predict.et_model <- function(model, newdata, env_noise = NULL,
                                  n_draws = 2000L,
                                  ci_levels = c(0.5, 0.8, 0.9, 0.95),
                                  n_perturb = NULL, ...) {

  fit         <- model$fit
  pred_names  <- .brms_pred_names(fit)
  n_perturb   <- if (is.null(n_perturb)) min(500L, n_draws) else as.integer(n_perturb)

  # Resolve environmental noise SDs
  noise_sds <- .resolve_env_noise(env_noise, pred_names, newdata)

  # --- 1. Posterior predictive (full uncertainty) ---
  pp <- brms::posterior_predict(fit, newdata = newdata, ndraws = n_draws)

  # --- 2. Linear predictor (parameter uncertainty only) ---
  lp <- brms::posterior_linpred(fit, newdata = newdata, ndraws = n_draws)

  # --- 3. Sigma draws ---
  draws_mat   <- .brms_draws_matrix(fit, max_draws = max(n_draws, n_perturb))
  sigma_draws <- if ("sigma" %in% colnames(draws_mat)) {
    draws_mat[seq_len(n_draws), "sigma"]
  } else {
    rep(NA_real_, n_draws)
  }

  # --- 4. Perturbed linear predictor (environmental uncertainty) ---
  # Short-circuit: if no predictor has noise, lp_perturbed == lp exactly,
  # so env_var will be zero without running the perturbation loop.
  lp_perturbed <- if (all(noise_sds == 0, na.rm = TRUE)) {
    lp[seq_len(n_perturb), , drop = FALSE]
  } else {
    .compute_lp_perturbed(
      draws_mat  = draws_mat[seq_len(n_perturb), , drop = FALSE],
      newdata    = newdata,
      pred_names = pred_names,
      noise_sds  = noise_sds
    )
  }

  # --- 5. Credible intervals ---
  ci_df <- .compute_ci(pp, ci_levels)

  # --- 6. Decomposition ---
  decomp <- .decompose_from_arrays(
    pp           = pp,
    lp           = lp,
    lp_perturbed = lp_perturbed,
    sigma_draws  = sigma_draws
  )

  structure(
    list(
      posterior_predict = pp,
      posterior_linpred = lp,
      lp_perturbed      = lp_perturbed,
      sigma_draws       = sigma_draws,
      credible_intervals = ci_df,
      decomposition     = decomp,
      newdata           = newdata,
      model             = model,
      env_noise         = env_noise,
      n_draws           = n_draws
    ),
    class = "et_prediction"
  )
}

# ============================================================
# Method: et_model_list (grouped models)
# ============================================================

#' @export
et_predict.et_model_list <- function(model, newdata, env_noise = NULL,
                                      n_draws = 2000L,
                                      ci_levels = c(0.5, 0.8, 0.9, 0.95),
                                      n_perturb = NULL, ...) {

  grouping <- model$grouping
  if (!grouping %in% colnames(newdata)) {
    stop("grouping column '", grouping, "' not found in newdata.")
  }

  groups <- names(model$models)
  preds  <- vector("list", length(groups))
  names(preds) <- groups

  for (g in groups) {
    m <- model$models[[g]]
    if (is.null(m)) {
      .et_warn("Skipping group ", g, " (model is NULL)")
      next
    }
    sub_nd <- newdata[newdata[[grouping]] == g, , drop = FALSE]
    if (nrow(sub_nd) == 0) {
      .et_warn("No newdata rows for group ", g, " — skipping")
      next
    }

    preds[[g]] <- tryCatch(
      et_predict.et_model(
        model     = m,
        newdata   = sub_nd,
        env_noise = env_noise,
        n_draws   = n_draws,
        ci_levels = ci_levels,
        n_perturb = n_perturb,
        ...
      ),
      error = function(e) {
        .et_error("Prediction failed for group ", g, ": ", e$message)
        NULL
      }
    )
  }

  structure(
    list(
      predictions = preds,
      grouping    = grouping,
      newdata     = newdata
    ),
    class = "et_prediction_list"
  )
}

# ============================================================
# Internal helpers
# ============================================================

# Compute lp for environmentally perturbed predictor values.
# draws_mat : [n_perturb x n_params] posterior draws
# newdata   : data.frame of predictors
# pred_names: character vector
# noise_sds : named numeric vector of per-predictor noise SDs
.compute_lp_perturbed <- function(draws_mat, newdata, pred_names, noise_sds) {
  n_perturb <- nrow(draws_mat)
  n_obs     <- nrow(newdata)
  lp_mat    <- matrix(NA_real_, n_perturb, n_obs)

  intercept_col <- "b_Intercept"
  beta_cols     <- paste0("b_", pred_names)
  avail_betas   <- intersect(beta_cols, colnames(draws_mat))
  avail_preds   <- sub("^b_", "", avail_betas)

  if (length(avail_betas) == 0) {
    .et_warn("No beta columns found in draws matrix — returning zero LP.")
    return(matrix(0, n_perturb, n_obs))
  }

  int_vals <- if (intercept_col %in% colnames(draws_mat)) {
    draws_mat[, intercept_col]
  } else {
    rep(0, n_perturb)
  }

  for (i in seq_len(n_perturb)) {
    nd_p <- newdata
    for (p in pred_names) {
      sd_p <- noise_sds[p]
      if (!is.na(sd_p) && sd_p > 0) {
        nd_p[[p]] <- nd_p[[p]] + stats::rnorm(n_obs, 0, sd_p)
      }
    }
    xmat   <- as.matrix(nd_p[, avail_preds, drop = FALSE])
    betas  <- draws_mat[i, avail_betas]
    lp_mat[i, ] <- int_vals[i] + xmat %*% betas
  }

  lp_mat
}

# Build credible interval data.frame from posterior predictive matrix.
# pp: [n_draws x n_obs]
.compute_ci <- function(pp, ci_levels) {
  n_obs <- ncol(pp)
  rows  <- vector("list", length(ci_levels))

  for (k in seq_along(ci_levels)) {
    lv    <- ci_levels[k]
    alpha <- 1 - lv
    lower  <- apply(pp, 2, stats::quantile, probs = alpha / 2)
    upper  <- apply(pp, 2, stats::quantile, probs = 1 - alpha / 2)
    med    <- apply(pp, 2, stats::median)
    rows[[k]] <- data.frame(
      row_id   = seq_len(n_obs),
      ci_level = lv,
      lower    = lower,
      median   = med,
      upper    = upper,
      width    = upper - lower,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

# Core uncertainty decomposition from arrays.
# Returns a data.frame with one row per observation.
.decompose_from_arrays <- function(pp, lp, lp_perturbed, sigma_draws) {
  n_obs <- ncol(pp)
  data.frame(
    obs_id       = seq_len(n_obs),
    total_var    = apply(pp, 2, stats::var),
    param_var    = apply(lp, 2, stats::var),
    env_var      = pmax(
      apply(lp_perturbed, 2, stats::var, na.rm = TRUE) -
        apply(lp[seq_len(nrow(lp_perturbed)), , drop = FALSE], 2, stats::var),
      0
    ),
    residual_var = mean(sigma_draws^2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# S3 methods for et_prediction
# ============================================================

#' @export
print.et_prediction <- function(x, ...) {
  cat("ErrorTracer prediction (et_prediction)\n")
  cat("  Observations :", ncol(x$posterior_predict), "\n")
  cat("  Draws        :", nrow(x$posterior_predict), "\n")
  cat("  CI levels    :", paste(unique(x$credible_intervals$ci_level), collapse = ", "), "\n")
  decomp <- x$decomposition
  cat("  Mean var decomposition (across observations):\n")
  cat(sprintf("    Parameter  : %.4f\n", mean(decomp$param_var)))
  cat(sprintf("    Environmental: %.4f\n", mean(decomp$env_var)))
  cat(sprintf("    Residual   : %.4f\n", mean(decomp$residual_var)))
  cat(sprintf("    Total      : %.4f\n", mean(decomp$total_var)))
  invisible(x)
}

#' @export
print.et_prediction_list <- function(x, ...) {
  cat("ErrorTracer grouped predictions (et_prediction_list)\n")
  cat("  Grouping :", x$grouping, "\n")
  n_ok <- sum(!vapply(x$predictions, is.null, logical(1)))
  cat("  Groups   :", n_ok, "/", length(x$predictions), "predicted\n")
  invisible(x)
}
