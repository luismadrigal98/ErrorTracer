# R/predict.R -- et_predict(): posterior prediction with uncertainty propagation

# ******************************************************************************
# S3 generic
# ______________________________________________________________________________

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
#'     \item \code{NULL} (default): no environmental noise.
#'     \item A single numeric: applied as a fraction of each predictor's
#'       empirical SD in \code{newdata} (e.g.\ \code{0.1} means 10\% noise,
#'       constant across all observations).
#'     \item A named \code{list} or named numeric vector with one scalar per
#'       predictor: constant absolute noise SD per predictor, e.g.\
#'       \code{list(Tmean = 0.5, PPT = 10)}.
#'     \item A named \code{list} where each entry is a \strong{numeric vector
#'       of length \code{nrow(newdata)}}: \emph{time-varying} (per-row) noise
#'       SDs.  Use this when predictor uncertainty grows with forecast horizon,
#'       e.g.\ from a GCM ensemble spread that increases over time:
#'       \code{list(Tmean = 0.30 + 0.01 * (years - base_year))}.
#'       Entries not supplied default to zero (no noise for that predictor).
#'   }
#' @param env_cov Correlation structure of the environmental noise.  The
#'   \emph{magnitudes} of the noise come from \code{env_noise}; \code{env_cov}
#'   supplies the \emph{correlation} between predictors, so that a perturbation
#'   on row \eqn{i} is drawn from
#'   \eqn{\mathcal{N}(0,\; D_i R D_i)} with
#'   \eqn{D_i = \mathrm{diag}(\sigma_{i1}, \ldots, \sigma_{ip})} and
#'   \eqn{R} = the correlation matrix.  One of:
#'   \itemize{
#'     \item \code{NULL} (default): independent noise, \eqn{R = I} --- equivalent
#'       to ErrorTracer behaviour prior to this feature and the right choice
#'       when predictor measurement errors are genuinely independent (e.g.\
#'       separate instruments on unrelated variables).
#'     \item \code{"empirical"}: compute the correlation of the predictor
#'       columns in the \strong{training data} (\code{model$data}).  Use this
#'       when predictor \emph{errors} are expected to inherit the correlation
#'       structure of the predictors themselves --- e.g.\ temperature and
#'       humidity that co-vary in the underlying climate system.
#'     \item \code{"newdata"}: compute the correlation of the predictor columns
#'       in \code{newdata}.  Useful when the forecast window has a different
#'       covariance structure than training (e.g.\ scenario runs).
#'     \item A numeric \eqn{p \times p} matrix with \code{dimnames} matching the
#'       model's predictors.  Entries with an off-diagonal exceeding 1 are
#'       rescaled to a correlation matrix.  Use this to supply an independent
#'       estimate of the \emph{error} correlation structure (e.g.\ from a
#'       reanalysis product or a sensor covariance report).
#'   }
#'   A correlation derived from training data is a working assumption: the
#'   structure of the \emph{errors} is assumed to mirror the structure of the
#'   \emph{values}.  When this is implausible, pass a matrix directly.
#' @param env_dist Distributional form of the per-predictor noise.  The
#'   \code{env_noise} SDs set the \emph{magnitude} of the perturbation;
#'   \code{env_dist} sets its \emph{shape}.  For every distribution other than
#'   \code{"gaussian"}, the noise is calibrated so that (approximately)
#'   \eqn{E[\tilde x] = x} and \eqn{\mathrm{Var}[\tilde x] = \sigma^2}, using a
#'   Gaussian copula to honour \code{env_cov}.  One of:
#'   \itemize{
#'     \item \code{NULL} (default): \code{"gaussian"} for every predictor ---
#'       additive Gaussian noise, legacy behaviour.
#'     \item A single string (\code{"gaussian"}, \code{"lognormal"},
#'       \code{"gamma"}, \code{"beta"}): applied to all predictors.
#'     \item A named list / character vector with one entry per predictor to
#'       override the default, e.g.\
#'       \code{list(PPT = "gamma", tmax = "gaussian")}.
#'   }
#'   Distributions:
#'   \describe{
#'     \item{\code{"gaussian"}}{Additive normal noise (\eqn{\tilde x = x + \varepsilon},
#'       \eqn{\varepsilon \sim N(0, \sigma^2)}).  Appropriate for symmetric
#'       measurement error on a continuous, potentially negative scale
#'       (temperature, anomalies).}
#'     \item{\code{"lognormal"}}{Multiplicative noise: \eqn{\log \tilde x \sim
#'       N(\log x - s^2/2,\; s^2)} with \eqn{s^2 = \log(1 + (\sigma/x)^2)}.
#'       Preserves positivity; right-tail skewed.  Natural for strictly
#'       positive continuous variables whose error scales with magnitude
#'       (e.g.\ enzyme activity, biomass).  Rows with \eqn{x \le 0} are
#'       left unperturbed.}
#'     \item{\code{"gamma"}}{\eqn{\tilde x \sim \mathrm{Gamma}(\mathrm{shape} =
#'       (x/\sigma)^2,\; \mathrm{rate} = x/\sigma^2)}.  Positive support,
#'       right-skewed, analytic mean/variance match.  Natural for precipitation,
#'       rates, and other non-negative continuous variables.  Rows with
#'       \eqn{x \le 0} are left unperturbed.}
#'     \item{\code{"beta"}}{\eqn{\tilde x \sim \mathrm{Beta}(\alpha,\beta)} with
#'       \eqn{\alpha + \beta = x(1-x)/\sigma^2 - 1}.  Support in \eqn{(0,1)};
#'       appropriate for proportions and probabilities (allele frequencies,
#'       presence rates).  Rows with \eqn{x \not\in (0,1)} or
#'       \eqn{\sigma^2 \ge x(1-x)} are left unperturbed.}
#'   }
#'   Correlation (\code{env_cov}) is applied to the latent standard-normal
#'   draws before the marginal quantile transform, so rank correlations are
#'   preserved across distributions.
#' @param n_draws Integer.  Number of posterior draws to use (default 2000;
#'   capped at the number of draws available in the fit).
#' @param ci_levels Numeric vector.  Credible interval levels to compute
#'   (default \code{c(0.5, 0.8, 0.9, 0.95)}).
#' @param n_perturb Integer.  Number of draws used for the environmental
#'   perturbation step (default \code{min(500, n_draws)}).  Reducing this
#'   speeds up computation.
#' @param interval_type Character.  Which draws to use when computing credible
#'   intervals:
#'   \itemize{
#'     \item \code{"predictive"} (default): draws from \code{posterior_predict},
#'       which include sigma (residual noise).  Use this when forecasting
#'       \strong{individual observations} — e.g. a single population's allele
#'       frequency, one site's ozone reading on a specific day.
#'     \item \code{"linpred"}: draws from \code{posterior_linpred}, which
#'       capture only parameter uncertainty (no sigma).  Use this when
#'       forecasting the \strong{mean response} — e.g. the expected ozone
#'       across many similar days, or mean Δf across replicate populations.
#'       These intervals are always narrower; they will under-cover individual
#'       observations unless sigma is negligible.
#'   }
#'   The decomposition components and \code{posterior_predict} /
#'   \code{posterior_linpred} matrices are always computed regardless of this
#'   setting.
#' @param include_env_in_ci Logical.  When \code{TRUE} and
#'   \code{interval_type = "predictive"}, credible intervals are constructed
#'   from \strong{environmentally inflated} draws
#'   \eqn{\tilde y = \tilde{\mathrm{lp}} + \varepsilon}, with
#'   \eqn{\tilde{\mathrm{lp}}} the perturbed linear predictor and
#'   \eqn{\varepsilon \sim N(0, \sigma^2)} using posterior draws of
#'   \eqn{\sigma}.  This folds the environmental uncertainty component back
#'   into the CI, which is typically what you want for sensitivity analyses
#'   or whenever the reported interval should cover predictor-measurement
#'   error.  When \code{FALSE} (default, backward compatible), CIs are based
#'   on \code{posterior_predict} only — parameter + residual, without
#'   predictor noise.
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
#'   \item{\code{env_cov}}{The \eqn{p \times p} correlation matrix actually
#'     used for perturbation (identity for \code{env_cov = NULL}).}
#'   \item{\code{env_dist}}{Named character vector mapping each predictor to
#'     the distribution actually used for its perturbation.}
#' }
#'
#' @seealso \code{\link{decompose_uncertainty}}, \code{\link{shelf_life}},
#'   \code{\link{et_calibrate}}
#' @export
et_predict <- function(model, newdata, env_noise = NULL, env_cov = NULL,
                        env_dist = NULL,
                        n_draws = 2000L, ci_levels = c(0.5, 0.8, 0.9, 0.95),
                        n_perturb = NULL,
                        interval_type = c("predictive", "linpred"),
                        include_env_in_ci = FALSE, ...) {
  UseMethod("et_predict")
}

# ******************************************************************************
# Method: et_model (single model)
# ______________________________________________________________________________

#' @export
et_predict.et_model <- function(model, newdata, env_noise = NULL,
                                  env_cov = NULL,
                                  env_dist = NULL,
                                  n_draws = 2000L,
                                  ci_levels = c(0.5, 0.8, 0.9, 0.95),
                                  n_perturb = NULL,
                                  interval_type = c("predictive", "linpred"),
                                  include_env_in_ci = FALSE,
                                  ...) {

  interval_type <- match.arg(interval_type)
  fit <- model$fit
  pred_names <- .brms_pred_names(fit)
  n_perturb <- if (is.null(n_perturb)) min(500L, n_draws) else as.integer(n_perturb)

  # For EIV-fit models, posterior_predict() on newdata requires each se_<pred>
  # column present. Copy them over from the training data (recycled to the
  # length of newdata) when the user hasn't supplied them.
  if (!is.null(model$eiv_spec)) {
    for (p in names(model$eiv_spec)) {
      se_col <- paste0("se_", p)
      if (!se_col %in% colnames(newdata)) {
        v <- as.numeric(model$eiv_spec[[p]])
        if (length(v) == 1L) v <- rep(v, nrow(newdata))
        if (length(v) != nrow(newdata)) {
          stop("eiv_spec[['", p, "']] has length ", length(v),
               " but newdata has ", nrow(newdata), " row(s). ",
               "Supply a se_", p, " column in newdata explicitly.")
        }
        newdata[[se_col]] <- v
      }
    }
  }

  # Resolve environmental noise SDs, correlation structure, and per-predictor
  # perturbation distribution.
  noise_sds <- .resolve_env_noise(env_noise, pred_names, newdata)
  cor_mat   <- .resolve_env_cor(env_cov, pred_names,
                                training_data = model$data,
                                newdata       = newdata)
  dist_spec <- .resolve_env_dist(env_dist, pred_names)

  # If the model was fit with errors-in-variables (me() terms), the fit's
  # beta posteriors already absorb predictor measurement noise; perturbing
  # those predictors again would double-count. Warn, so the user can
  # null-out env_noise for EIV-modelled predictors.
  eiv_preds <- names(model$eiv_spec %||% list())
  if (length(eiv_preds) && !is.null(env_noise)) {
    overlap <- intersect(eiv_preds, names(noise_sds))
    overlap <- overlap[vapply(overlap, function(p) any(noise_sds[[p]] != 0),
                              logical(1))]
    if (length(overlap)) {
      .et_warn("env_noise is non-zero for predictor(s) also modelled ",
               "under eiv (", paste(overlap, collapse = ", "),
               "); env_var may double-count predictor uncertainty.")
    }
  }

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
  # Short-circuit: if every predictor's noise is identically zero for every
  # observation, lp_perturbed == lp and env_var will be zero -- skip the loop.
  all_zero_noise <- all(vapply(noise_sds,
                               function(v) all(v == 0, na.rm = TRUE),
                               logical(1)))
  lp_perturbed <- if (all_zero_noise) {
    lp[seq_len(n_perturb), , drop = FALSE]
  } else {
    .compute_lp_perturbed(
      draws_mat  = draws_mat[seq_len(n_perturb), , drop = FALSE],
      newdata    = newdata,
      pred_names = pred_names,
      noise_sds  = noise_sds,   # named list of per-obs vectors
      cor_mat    = cor_mat,     # p x p correlation over pred_names
      dist_spec  = dist_spec    # named character: predictor -> distribution
    )
  }

  # --- 5. Credible intervals ---
  # Route to posterior_predict (includes sigma) or posterior_linpred (no sigma)
  # depending on interval_type. The underlying matrices are always stored for
  # decomposition regardless.
  # When include_env_in_ci = TRUE and interval_type is predictive, rebuild
  # predictive draws as lp_perturbed + sigma * N(0,1) so predictor noise is
  # folded back into the credible interval.
  ci_draws <- if (interval_type == "predictive") {
    if (isTRUE(include_env_in_ci) && !all_zero_noise) {
      .inflate_env_predictive(lp_perturbed, sigma_draws)
    } else {
      pp
    }
  } else {
    lp
  }
  ci_df <- .compute_ci(ci_draws, ci_levels)

  # --- 6. Decomposition ---
  decomp <- .decompose_from_arrays(
    pp           = pp,
    lp           = lp,
    lp_perturbed = lp_perturbed,
    sigma_draws  = sigma_draws
  )

  structure(
    list(
      posterior_predict  = pp,
      posterior_linpred  = lp,
      lp_perturbed       = lp_perturbed,
      sigma_draws        = sigma_draws,
      credible_intervals = ci_df,
      decomposition      = decomp,
      newdata            = newdata,
      model              = model,
      env_noise          = env_noise,
      env_cov            = cor_mat,
      env_dist           = dist_spec,
      n_draws            = n_draws,
      interval_type      = interval_type,
      include_env_in_ci  = isTRUE(include_env_in_ci)
    ),
    class = "et_prediction"
  )
}

# ******************************************************************************
# Method: et_model_list (grouped models)
# ______________________________________________________________________________

#' @export
et_predict.et_model_list <- function(model, newdata, env_noise = NULL,
                                      env_cov = NULL,
                                      env_dist = NULL,
                                      n_draws = 2000L,
                                      ci_levels = c(0.5, 0.8, 0.9, 0.95),
                                      n_perturb = NULL,
                                      interval_type = c("predictive", "linpred"),
                                      include_env_in_ci = FALSE,
                                      ...) {
  interval_type <- match.arg(interval_type)

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
      .et_warn("No newdata rows for group ", g, " -- skipping")
      next
    }

    preds[[g]] <- tryCatch(
      et_predict.et_model(
        model             = m,
        newdata           = sub_nd,
        env_noise         = env_noise,
        env_cov           = env_cov,
        env_dist          = env_dist,
        n_draws           = n_draws,
        ci_levels         = ci_levels,
        n_perturb         = n_perturb,
        interval_type     = interval_type,
        include_env_in_ci = include_env_in_ci,
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

# ******************************************************************************
# Internal helpers
# ______________________________________________________________________________

# Compute lp for environmentally perturbed predictor values.
# draws_mat : [n_perturb x n_params] posterior draws
# newdata   : data.frame of predictors
# pred_names: character vector
# noise_sds : named list; each element is a numeric vector of length n_obs
#             giving the per-observation noise SD for that predictor.
# cor_mat   : p x p correlation matrix over pred_names (identity = independent).
#             Used to draw correlated standard-normal latents; marginal
#             distributions are then applied per predictor via dist_spec.
# dist_spec : named character vector (predictor -> distribution). Controls
#             the marginal perturbation form (gaussian / lognormal / gamma
#             / beta) via .perturb_predictor(). The correlation in cor_mat is
#             applied on the Gaussian copula so rank correlations are
#             preserved across distributions.
.compute_lp_perturbed <- function(draws_mat, newdata, pred_names, noise_sds,
                                   cor_mat = NULL, dist_spec = NULL) {
  n_perturb <- nrow(draws_mat)
  n_obs     <- nrow(newdata)
  lp_mat    <- matrix(NA_real_, n_perturb, n_obs)

  intercept_col <- "b_Intercept"
  beta_cols     <- paste0("b_", pred_names)
  avail_betas   <- intersect(beta_cols, colnames(draws_mat))
  avail_preds   <- sub("^b_", "", avail_betas)

  if (length(avail_betas) == 0) {
    .et_warn("No beta columns found in draws matrix -- returning zero LP.")
    return(matrix(0, n_perturb, n_obs))
  }

  int_vals <- if (intercept_col %in% colnames(draws_mat)) {
    draws_mat[, intercept_col]
  } else {
    rep(0, n_perturb)
  }

  # Cholesky of the correlation sub-matrix restricted to available predictors.
  # L is upper triangular with t(L) %*% L == R, so rows of Z %*% L with
  # Z ~ N(0, I) have covariance R.
  L <- NULL
  if (!is.null(cor_mat)) {
    R_sub <- cor_mat[avail_preds, avail_preds, drop = FALSE]
    if (!isTRUE(all.equal(R_sub, diag(length(avail_preds)),
                          check.attributes = FALSE))) {
      L <- tryCatch(chol(R_sub),
                    error = function(e) {
                      .et_warn("Cholesky of env_cov failed (", e$message,
                               "); falling back to independent noise.")
                      NULL
                    })
    }
  }

  # Stack per-predictor per-obs SDs into an n_obs x p matrix for elementwise
  # scaling of the standardized draws.
  sd_mat <- do.call(cbind, lapply(avail_preds, function(p) {
    v <- noise_sds[[p]]
    if (is.null(v)) rep(0, n_obs) else v
  }))
  colnames(sd_mat) <- avail_preds

  # Fill in any missing distribution entries with "gaussian".
  if (is.null(dist_spec)) {
    dist_spec <- stats::setNames(rep("gaussian", length(avail_preds)),
                                 avail_preds)
  } else {
    missing_preds <- setdiff(avail_preds, names(dist_spec))
    if (length(missing_preds)) {
      dist_spec[missing_preds] <- "gaussian"
    }
  }

  # Original predictor values as a matrix for distributional perturbation.
  x_mat <- as.matrix(newdata[, avail_preds, drop = FALSE])

  for (i in seq_len(n_perturb)) {
    Z <- matrix(stats::rnorm(n_obs * length(avail_preds)),
                n_obs, length(avail_preds))
    if (!is.null(L)) Z <- Z %*% L

    xmat_p <- x_mat
    for (k in seq_along(avail_preds)) {
      p <- avail_preds[k]
      xmat_p[, k] <- .perturb_predictor(
        x     = x_mat[, k],
        sigma = sd_mat[, k],
        Z_std = Z[, k],
        dist  = dist_spec[[p]]
      )
    }

    betas       <- draws_mat[i, avail_betas]
    lp_mat[i, ] <- int_vals[i] + xmat_p %*% betas
  }

  lp_mat
}

# Fold environmental noise into the predictive draws by adding fresh
# residual noise on top of the perturbed linear predictor. The result is a
# matrix with the same shape as lp_perturbed whose column variance equals
# approximately var(lp_perturbed) + mean(sigma_draws^2), so CIs built from
# it cover parameter + env + residual uncertainty.
.inflate_env_predictive <- function(lp_perturbed, sigma_draws) {
  n_p  <- nrow(lp_perturbed)
  n_o  <- ncol(lp_perturbed)
  s    <- sigma_draws[seq_len(n_p)]
  s[is.na(s)] <- 0
  # rnorm with vector sd is recycled element-wise; wrapping in matrix(, nrow=n_p)
  # then aligns row i with sigma i since R fills column-major.
  eps  <- matrix(stats::rnorm(n_p * n_o, sd = rep(s, n_o)), nrow = n_p)
  lp_perturbed + eps
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

# ******************************************************************************
# S3 methods for et_prediction
# ______________________________________________________________________________

#' @export
print.et_prediction <- function(x, ...) {
  cat("ErrorTracer prediction (et_prediction)\n")
  cat("  Observations  :", ncol(x$posterior_predict), "\n")
  cat("  Draws         :", nrow(x$posterior_predict), "\n")
  cat("  CI levels     :", paste(unique(x$credible_intervals$ci_level), collapse = ", "), "\n")
  cat("  Interval type :", if (is.null(x$interval_type)) "predictive" else x$interval_type, "\n")
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
