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
#' @param env_ensemble Optional pathwise driver ensemble for forecasting under
#'   uncertain future covariates.  A \strong{named list}, one entry per uncertain
#'   driver predictor, where each entry is an \eqn{M \times H} matrix of
#'   \eqn{M} scenario trajectories over the \eqn{H = \mathrm{nrow(newdata)}}
#'   forecast steps (rows are scenarios, columns are horizon steps).  Unlike
#'   \code{env_noise} --- which draws \emph{independent} per-observation jitter
#'   that averages out and leaves the driver variance flat with lead time ---
#'   each posterior draw here is paired with a whole \emph{trajectory}, so the
#'   driver uncertainty is temporally coherent and \strong{accumulates with lead
#'   time} as a real driver ensemble (e.g.\ CMIP members, or resampled
#'   climatology) fans out.  Predictors absent from \code{env_ensemble} are held
#'   at their \code{newdata} values.  When supplied, \code{env_ensemble}
#'   supersedes \code{env_noise} for the perturbation and is folded into the
#'   credible interval by default.  This is the recommended way to build a
#'   genuine \emph{reforecast} in which only information available at the issue
#'   date is used.
#' @param n_draws Integer.  Number of posterior draws to use (default 2000;
#'   capped at the number of draws available in the fit).
#' @param ci_levels Numeric vector.  Credible interval levels to compute
#'   (default \code{c(0.5, 0.8, 0.9, 0.95)}).
#' @param n_perturb Integer.  Number of posterior draws used for the
#'   environmental perturbation step (default \code{min(500, n_draws)}).
#'   Reducing this speeds up computation.
#' @param n_env_draws Integer.  Number of independent environmental
#'   perturbations averaged \emph{per posterior draw} when estimating
#'   \code{env_var} (default \code{1}).  Increasing this reduces Monte Carlo
#'   noise on the environmental-variance estimate at the cost of proportional
#'   computation.  The decomposition always reports a Monte Carlo SE
#'   (\code{v_env_mcse}) alongside \code{env_var}; it decreases roughly like
#'   \eqn{1/\sqrt{n\_env\_draws \cdot n\_perturb}}.
#' @param interval_type Character.  Which draws to use when computing credible
#'   intervals:
#'   \itemize{
#'     \item \code{"predictive"} (default): draws from \code{posterior_predict},
#'       which include sigma (residual noise).  Use this when forecasting
#'       \strong{individual observations} --- e.g. a single population's allele
#'       frequency, one site's ozone reading on a specific day.
#'     \item \code{"linpred"}: draws from \code{posterior_linpred}, which
#'       capture only parameter uncertainty (no sigma).  Use this when
#'       forecasting the \strong{mean response} --- e.g. the expected ozone
#'       across many similar days, or mean delta f across replicate populations.
#'       These intervals are always narrower; they will under-cover individual
#'       observations unless sigma is negligible.
#'   }
#'   The decomposition components and \code{posterior_predict} /
#'   \code{posterior_linpred} matrices are always computed regardless of this
#'   setting.
#' @param include_env_in_ci Logical or \code{NULL}.  Controls whether
#'   environmental (predictor) uncertainty is folded into the credible
#'   interval when \code{interval_type = "predictive"}.  When \code{TRUE},
#'   intervals are built from \strong{environmentally inflated} draws in which
#'   each posterior-predictive residual is re-centred on the env-perturbed
#'   response-scale mean, \eqn{\tilde y = \tilde\mu + (y^{pp} - \mu)}
#'   (family-general; for Gaussian identity this equals
#'   \eqn{\tilde{\mathrm{lp}} + \varepsilon}, \eqn{\varepsilon\sim N(0,\sigma^2)}).
#'   The interval variance is then param + env + residual, coherent with the
#'   env-inclusive \code{total_var}.  The default \code{NULL} means
#'   \strong{auto}: env is folded in whenever the caller supplied non-zero
#'   \code{env_noise}, and omitted otherwise.  Pass \code{FALSE} to force a
#'   parameter + residual interval even when \code{env_noise} is given.
#' @param ... Passed to methods.
#'
#' @return An \code{et_prediction} object (list) containing:
#' \describe{
#'   \item{\code{posterior_predict}}{Matrix \code{[n_draws x n_obs]}: full
#'     posterior predictive draws (parameter + residual uncertainty).}
#'   \item{\code{posterior_linpred}}{Matrix \code{[n_draws x n_obs]}: linear
#'     predictor draws on the \strong{link scale} (parameter uncertainty only).}
#'   \item{\code{lp_perturbed}}{Matrix \code{[n_perturb x n_obs]}: linear
#'     predictor on the link scale computed on environmentally perturbed inputs.}
#'   \item{\code{sigma_draws}}{Numeric vector: posterior draws of sigma
#'     (\code{NA} for families without a sigma parameter, e.g.\ Binomial).}
#'   \item{\code{credible_intervals}}{data.frame with columns
#'     \code{row_id, ci_level, lower, median, upper, width}.}
#'   \item{\code{decomposition}}{data.frame with columns
#'     \code{obs_id, total_var, param_var, env_var, v_env_mcse, residual_var}.
#'     All components are on the \strong{response scale}.
#'     \code{v_env_mcse} is the Monte Carlo SE of \code{env_var}.
#'     \code{residual_var} is per-observation for non-Gaussian families
#'     (e.g.\ Binomial: \eqn{E_s[\mu^{(s)}(1-\mu^{(s)})]}) and constant
#'     for Gaussian.}
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
                        env_dist = NULL, env_ensemble = NULL,
                        n_draws = 2000L, ci_levels = c(0.5, 0.8, 0.9, 0.95),
                        n_perturb = NULL,
                        n_env_draws = 1L,
                        interval_type = c("predictive", "linpred"),
                        include_env_in_ci = NULL, ...) {
  UseMethod("et_predict")
}

# ******************************************************************************
# Method: et_model (single model)
# ______________________________________________________________________________

#' @export
et_predict.et_model <- function(model, newdata, env_noise = NULL,
                                  env_cov = NULL,
                                  env_dist = NULL,
                                  env_ensemble = NULL,
                                  n_draws = 2000L,
                                  ci_levels = c(0.5, 0.8, 0.9, 0.95),
                                  n_perturb = NULL,
                                  n_env_draws = 1L,
                                  interval_type = c("predictive", "linpred"),
                                  include_env_in_ci = NULL,
                                  ...) {

  interval_type <- match.arg(interval_type)
  n_env_draws   <- max(1L, as.integer(n_env_draws))
  fit           <- model$fit
  pred_names    <- .brms_pred_names(fit)
  n_perturb     <- if (is.null(n_perturb)) min(500L, n_draws) else as.integer(n_perturb)

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

  # Detect Gaussian identity link once; used to avoid a redundant
  # posterior_linpred(transform=TRUE) call and to stay backward-compatible.
  is_gauss_id <- .is_gaussian_identity(fit$family)

  # Detect a residual autocorrelation structure once. When present, the
  # posterior predictive must be generated by iterating the fitted process
  # forward across the forecast horizon (see .posterior_predict_forecast),
  # otherwise its variance does not accumulate with lead time.
  has_autocor <- .formula_has_autocor(fit$formula)

  # Cap the number of draws at what the fit actually holds, then use an
  # explicit, shared set of draw ids so posterior_predict and
  # posterior_linpred are aligned draw-for-draw (required for the residual
  # telescoping in .inflate_env_predictive and to avoid brms's random
  # sub-sampling returning mismatched draws across the two calls).
  n_avail <- brms::ndraws(fit)
  if (n_draws > n_avail) {
    .et_warn("n_draws (", n_draws, ") exceeds the ", n_avail,
             " draws in the fit; using ", n_avail, ".")
    n_draws <- n_avail
  }
  # The perturbation subsample cannot exceed the (capped) number of draws,
  # since mu_draws_sub indexes the first n_perturb rows of the draw matrices.
  n_perturb <- min(n_perturb, n_draws)
  draw_ids  <- seq_len(n_draws)

  # --- 1. Posterior predictive (full uncertainty, response scale) ---
  # For autocorrelation models, iterate the process across the horizon by
  # predicting on the training series continued with the forecast rows
  # (response = NA); brms then accumulates the AR/MA variance with lead time
  # and carries the initial-condition uncertainty from the forecast origin.
  pp <- if (has_autocor) {
    .posterior_predict_forecast(fit, model$data, newdata, draw_ids)
  } else {
    brms::posterior_predict(fit, newdata = newdata, draw_ids = draw_ids)
  }

  # --- 2. Linear predictor draws (link scale, parameter uncertainty only) ---
  # The mean structure carries no autocorrelation (ar()/ma() act on the
  # residual), so the linear predictor is computed directly on newdata.
  lp <- brms::posterior_linpred(fit, newdata = newdata, draw_ids = draw_ids)

  # --- 3. Response-scale parameter draws ---
  # For Gaussian identity g^{-1}(eta) = eta, so mu_draws == lp — reuse it.
  # For all other families, apply the inverse link to get the response scale.
  mu_draws <- if (is_gauss_id) {
    lp
  } else {
    brms::posterior_linpred(fit, newdata = newdata, draw_ids = draw_ids,
                            transform = TRUE)
  }

  # --- 4. Draws matrix and dispersion parameters ---
  draws_mat   <- .brms_draws_matrix(fit, max_draws = max(n_draws, n_perturb))
  sigma_draws <- if ("sigma" %in% colnames(draws_mat)) {
    draws_mat[seq_len(n_draws), "sigma"]
  } else {
    rep(NA_real_, n_draws)
  }
  disp_draws <- .extract_disp_draws(draws_mat, n_draws)

  # --- 5. Perturbed linear predictor (environmental uncertainty, link scale) ---
  # Short-circuit: if every predictor's noise is identically zero for every
  # observation, lp_perturbed == lp and env_var will be zero -- skip the loop.
  all_zero_noise <- all(vapply(noise_sds,
                               function(v) all(v == 0, na.rm = TRUE),
                               logical(1)))

  # Pathwise driver-ensemble mode (T2): when env_ensemble is supplied, each
  # posterior draw is paired with a whole covariate TRAJECTORY, so driver
  # uncertainty is temporally coherent and accumulates with lead time (unlike
  # iid per-observation env_noise). It supersedes iid noise for the perturbation.
  has_ensemble <- !is.null(env_ensemble)
  if (has_ensemble) {
    env_ensemble <- .validate_env_ensemble(env_ensemble, pred_names, newdata)
  }

  # There is environmental (driver) uncertainty if iid noise is non-zero OR a
  # driver ensemble was supplied.
  has_env_uncertainty <- (!all_zero_noise) || has_ensemble

  # Resolve the auto (NULL) default for include_env_in_ci: fold environmental
  # uncertainty into the interval whenever there is any. This keeps the reported
  # forecast interval (and hence shelf life + calibration coverage) coherent
  # with the env-inclusive total_var. Users can still override with TRUE/FALSE.
  if (is.null(include_env_in_ci)) include_env_in_ci <- has_env_uncertainty

  lp_perturbed <- if (has_ensemble) {
    n_scen       <- nrow(env_ensemble[[1]])
    scenario_ids <- rep_len(seq_len(n_scen), n_perturb)  # balanced, deterministic
    .compute_lp_ensemble(
      draws_mat    = draws_mat[seq_len(n_perturb), , drop = FALSE],
      newdata      = newdata,
      pred_names   = pred_names,
      env_ensemble = env_ensemble,
      scenario_ids = scenario_ids
    )
  } else if (all_zero_noise) {
    lp[seq_len(n_perturb), , drop = FALSE]
  } else {
    .compute_lp_perturbed(
      draws_mat   = draws_mat[seq_len(n_perturb), , drop = FALSE],
      newdata     = newdata,
      pred_names  = pred_names,
      noise_sds   = noise_sds,
      cor_mat     = cor_mat,
      dist_spec   = dist_spec,
      n_env_draws = n_env_draws
    )
  }

  # --- 6. Response-scale perturbed draws (for decomposition) ---
  mu_perturbed <- if (is_gauss_id) {
    lp_perturbed
  } else {
    .apply_inv_link(lp_perturbed, fit$family)
  }
  mu_draws_sub <- mu_draws[seq_len(n_perturb), , drop = FALSE]

  # --- 7. Credible intervals ---
  # Route to posterior_predict (includes sigma) or posterior_linpred (no sigma)
  # depending on interval_type. When include_env_in_ci = TRUE and interval_type
  # is predictive, re-centre each posterior-predictive residual on the
  # env-perturbed mean so predictor noise is folded into the interval
  # (family-general; see .inflate_env_predictive).
  # Full predictive draws implied by the reported forecast, with environmental
  # (driver) uncertainty folded in when applicable. Single source of truth for
  # the predictive interval AND for PIT / calibration, so all three agree.
  predictive_draws <- if (isTRUE(include_env_in_ci) && has_env_uncertainty) {
    .inflate_env_predictive(mu_perturbed, pp, mu_draws)
  } else {
    pp
  }
  ci_draws <- if (interval_type == "predictive") predictive_draws else lp
  ci_df <- .compute_ci(ci_draws, ci_levels)

  # --- 8. Decomposition (all components on response scale) ---
  # has_autocor was detected above (used to route the posterior predictive).
  decomp <- .decompose_from_arrays(
    pp           = pp,
    mu_draws     = mu_draws,
    mu_perturbed = mu_perturbed,
    mu_draws_sub = mu_draws_sub,
    family       = fit$family,
    disp_draws   = disp_draws,
    has_autocor  = has_autocor
  )

  structure(
    list(
      posterior_predict  = pp,
      posterior_linpred  = lp,
      lp_perturbed       = lp_perturbed,
      predictive_draws   = predictive_draws,
      sigma_draws        = sigma_draws,
      credible_intervals = ci_df,
      decomposition      = decomp,
      newdata            = newdata,
      model              = model,
      env_noise          = env_noise,
      env_cov            = cor_mat,
      env_dist           = dist_spec,
      env_ensemble       = if (has_ensemble)
                             list(n_scenarios = nrow(env_ensemble[[1]]),
                                  drivers     = names(env_ensemble))
                           else NULL,
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
                                      env_ensemble = NULL,
                                      n_draws = 2000L,
                                      ci_levels = c(0.5, 0.8, 0.9, 0.95),
                                      n_perturb = NULL,
                                      n_env_draws = 1L,
                                      interval_type = c("predictive", "linpred"),
                                      include_env_in_ci = NULL,
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
        env_ensemble      = env_ensemble,
        n_draws           = n_draws,
        ci_levels         = ci_levels,
        n_perturb         = n_perturb,
        n_env_draws       = n_env_draws,
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
# draws_mat   : [n_perturb x n_params] posterior draws
# newdata     : data.frame of predictors
# pred_names  : character vector
# noise_sds   : named list; each element is a numeric vector of length n_obs
#               giving the per-observation noise SD for that predictor.
# cor_mat     : p x p correlation matrix over pred_names.
# dist_spec   : named character vector (predictor -> distribution).
# n_env_draws : number of perturbations to average per posterior draw.
#               When > 1, the within-draw MC noise on the linear predictor is
#               reduced by averaging k independent perturbations before taking
#               the cross-draw variance.
.compute_lp_perturbed <- function(draws_mat, newdata, pred_names, noise_sds,
                                   cor_mat = NULL, dist_spec = NULL,
                                   n_env_draws = 1L) {
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

  # Stack per-predictor per-obs SDs into an n_obs x p matrix.
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
    if (length(missing_preds)) dist_spec[missing_preds] <- "gaussian"
  }

  x_mat <- as.matrix(newdata[, avail_preds, drop = FALSE])
  p     <- length(avail_preds)

  for (i in seq_len(n_perturb)) {
    betas <- draws_mat[i, avail_betas]

    if (n_env_draws <= 1L) {
      # Single perturbation (original fast path)
      Z <- matrix(stats::rnorm(n_obs * p), n_obs, p)
      if (!is.null(L)) Z <- Z %*% L
      xmat_p <- x_mat
      for (k in seq_along(avail_preds)) {
        xmat_p[, k] <- .perturb_predictor(
          x     = x_mat[, k],
          sigma = sd_mat[, k],
          Z_std = Z[, k],
          dist  = dist_spec[[avail_preds[k]]]
        )
      }
      lp_mat[i, ] <- int_vals[i] + as.numeric(xmat_p %*% betas)
    } else {
      # Average n_env_draws perturbations; add the intercept once at the end.
      lp_acc <- rep(0, n_obs)
      for (j in seq_len(n_env_draws)) {
        Z <- matrix(stats::rnorm(n_obs * p), n_obs, p)
        if (!is.null(L)) Z <- Z %*% L
        xmat_p <- x_mat
        for (k in seq_along(avail_preds)) {
          xmat_p[, k] <- .perturb_predictor(
            x     = x_mat[, k],
            sigma = sd_mat[, k],
            Z_std = Z[, k],
            dist  = dist_spec[[avail_preds[k]]]
          )
        }
        lp_acc <- lp_acc + as.numeric(xmat_p %*% betas)
      }
      lp_mat[i, ] <- int_vals[i] + lp_acc / n_env_draws
    }
  }

  lp_mat
}

# Validate and normalise a driver ensemble (T2). Returns a named list of
# numeric [M x H] matrices (M scenarios x H = nrow(newdata) horizon steps),
# restricted to predictors that are actually in the model.
.validate_env_ensemble <- function(env_ensemble, pred_names, newdata) {
  if (!is.list(env_ensemble) || is.null(names(env_ensemble)) ||
      any(names(env_ensemble) == "")) {
    stop("env_ensemble must be a NAMED list of [M x H] matrices, one per ",
         "driver predictor (M scenarios by H = nrow(newdata) horizon steps).")
  }
  H <- nrow(newdata)

  unknown <- setdiff(names(env_ensemble), pred_names)
  if (length(unknown)) {
    .et_warn("env_ensemble contains predictor(s) not in model: ",
             paste(unknown, collapse = ", "), " -- ignored.")
    env_ensemble <- env_ensemble[intersect(names(env_ensemble), pred_names)]
  }
  if (length(env_ensemble) == 0L) {
    stop("env_ensemble has no drivers matching the model predictors.")
  }

  mats <- lapply(env_ensemble, function(e) {
    m <- as.matrix(e)
    storage.mode(m) <- "double"
    m
  })
  Ms <- vapply(mats, nrow, integer(1))
  Hs <- vapply(mats, ncol, integer(1))
  if (any(Hs != H)) {
    stop("each env_ensemble[[p]] must have ncol == nrow(newdata) (= ", H,
         "); a driver had ", paste(unique(Hs), collapse = ", "), " columns. ",
         "Rows are scenarios, columns are forecast/horizon steps.")
  }
  if (length(unique(Ms)) != 1L) {
    stop("all env_ensemble drivers must share the same number of scenarios ",
         "(rows); got ", paste(unique(Ms), collapse = ", "), ".")
  }
  if (unique(Ms) < 2L) {
    stop("env_ensemble needs at least 2 scenarios to carry driver uncertainty.")
  }
  mats
}

# Pathwise driver-ensemble linear predictor (T2). Instead of iid per-observation
# noise, each posterior draw is paired with a whole covariate TRAJECTORY drawn
# from a user-supplied driver ensemble (e.g. CMIP members or resampled
# climatology). Because a trajectory is temporally coherent -- a warm scenario
# stays warm across the whole horizon -- and ensembles fan out with lead time,
# the resulting driver variance ACCUMULATES with lead time, unlike iid jitter
# (which averages out and stays flat). Non-ensemble predictors are held at their
# newdata values.
#
# draws_mat    : [n_perturb x n_params] posterior draws (betas + intercept)
# newdata      : [H x .] forecast rows; supplies the fixed predictors
# pred_names   : model predictor names
# env_ensemble : named list; env_ensemble[[p]] is an [M x H] matrix of scenario
#                trajectories for driver predictor p
# scenario_ids : integer vector length n_perturb, the scenario paired with each
#                draw
.compute_lp_ensemble <- function(draws_mat, newdata, pred_names,
                                 env_ensemble, scenario_ids) {
  n_perturb <- nrow(draws_mat)
  H         <- nrow(newdata)

  beta_cols   <- paste0("b_", pred_names)
  avail_betas <- intersect(beta_cols, colnames(draws_mat))
  avail_preds <- sub("^b_", "", avail_betas)
  if (length(avail_betas) == 0L) {
    .et_warn("No beta columns found in draws matrix -- returning zero LP.")
    return(matrix(0, n_perturb, H))
  }

  int_vals <- if ("b_Intercept" %in% colnames(draws_mat)) {
    draws_mat[, "b_Intercept"]
  } else {
    rep(0, n_perturb)
  }

  X_base <- as.matrix(newdata[, avail_preds, drop = FALSE])  # H x P
  driver_names <- intersect(names(env_ensemble), avail_preds)

  lp_mat <- matrix(NA_real_, n_perturb, H)
  for (i in seq_len(n_perturb)) {
    Xi <- X_base
    m  <- scenario_ids[i]
    for (p in driver_names) {
      Xi[, p] <- env_ensemble[[p]][m, ]   # scenario m's trajectory for driver p
    }
    lp_mat[i, ] <- int_vals[i] + as.numeric(Xi %*% draws_mat[i, avail_betas])
  }
  lp_mat
}

# Forecast-aware posterior predictive for autocorrelation models
# (ar()/ma()/arma()/cosy()/unstr()/sar()/car()). brms accumulates the residual
# autocorrelation variance across a horizon only when the forecast rows
# CONTINUE the training series; predicting on the forecast rows alone yields a
# flat, non-accumulating predictive (the variance never grows with lead time).
# We therefore append the forecast rows (response = NA) to the training data,
# predict on the combined series, and return only the forecast columns. brms
# then iterates the fitted process forward — integrating over the posterior of
# the AR coefficient (so draws with |phi| >= 1 contribute super-linear growth)
# and carrying the initial-condition uncertainty from the forecast origin.
#
# fit           : brmsfit
# training_data : the data the model was fit on (model$data)
# newdata       : forecast rows (predictors + the ar() time index; the response
#                 column may be absent or NA)
# draw_ids      : integer draw indices (aligned with the other draw matrices)
.posterior_predict_forecast <- function(fit, training_data, newdata, draw_ids) {
  resp <- tryCatch(all.vars(brms::brmsterms(fit$formula)$respform),
                   error = function(e) NULL)
  if (is.null(resp) || length(resp) == 0L) {
    # Cannot identify the response: fall back to the plain (non-accumulating)
    # call rather than fail. Variance growth with lead time will be lost.
    .et_warn("Could not identify the response for an autocorrelation model; ",
             "forecast variance may not accumulate with lead time.")
    return(brms::posterior_predict(fit, newdata = newdata, draw_ids = draw_ids))
  }

  nd <- newdata
  # ALWAYS blank the response on the forecast rows: they are forecasts, so
  # their responses are unknown at the issue date. If newdata happens to carry
  # observed values (e.g. a hindcast for validation), letting brms see them
  # would condition on the future and yield smoothed, non-accumulating
  # intervals instead of a genuine forecast.
  for (r in resp) nd[[r]] <- NA
  # Align to the training columns; fill any still-missing columns with NA so
  # rbind() succeeds (a genuinely missing predictor surfaces downstream).
  miss <- setdiff(names(training_data), names(nd))
  for (m in miss) nd[[m]] <- NA
  nd <- nd[, names(training_data), drop = FALSE]

  combined <- rbind(training_data, nd)
  n_new    <- nrow(newdata)
  n_all    <- nrow(combined)

  pp_all <- brms::posterior_predict(fit, newdata = combined, draw_ids = draw_ids)
  pp_all[, (n_all - n_new + 1L):n_all, drop = FALSE]
}

# Fold environmental (predictor) uncertainty into the predictive draws by
# re-centring each posterior-predictive residual on the env-perturbed mean.
# Family-general: for draw s and observation i,
#   y_env^(s) = mu_perturbed^(s) + (pp^(s) - mu^(s))
# where (pp - mu) is the family residual implied by the UNPERTURBED posterior
# predictive (pp) and its response-scale mean (mu = mu_draws). This works for
# any brms family without a family-specific sampler: for Gaussian identity it
# reduces to mu_perturbed + sigma * N(0,1). The result has variance
# Var(mu_perturbed) + Var(residual) = param + env + residual, matching the
# variance decomposition, so the interval and total_var tell a consistent
# story. All arguments are on the RESPONSE scale.
#
# mu_perturbed : [n_perturb x n_obs] env-perturbed response-scale means
# pp           : [n_draws  x n_obs] unperturbed posterior predictive draws
# mu_draws     : [n_draws  x n_obs] unperturbed response-scale means
.inflate_env_predictive <- function(mu_perturbed, pp, mu_draws) {
  n_draws <- nrow(pp)
  n_p     <- nrow(mu_perturbed)
  # Recycle the env perturbation across ALL posterior draws so the returned
  # predictive draws cover the same n_draws as pp (and hence the same draws the
  # variance decomposition uses). Env perturbations are iid and exchangeable, so
  # recycling the pattern is valid, and the first n_p rows are unchanged.
  # Without this the interval / PIT / calibration would use only the first
  # n_perturb draws while the decomposition uses all n_draws -- an inconsistency
  # that is acute for autocorrelation forecasts whose far-lead draws are
  # heavy-tailed (a few near-non-stationary draws dominate the full sample but
  # not the truncated one).
  idx   <- rep_len(seq_len(n_p), n_draws)
  resid <- pp - mu_draws
  mu_perturbed[idx, , drop = FALSE] + resid
}

# Tail-robust (winsorized) column variance. Each column is clipped to its
# [trim, 1 - trim] quantiles before var() is taken, so a handful of extreme
# draws cannot dominate the estimate. This matters for the posterior-predictive
# variance of an autocorrelation forecast: when the AR coefficient's posterior
# has mass at |phi| >= 1, the k-step-ahead forecast variance is heavy-tailed and
# a single non-stationary draw can inflate the raw variance by an order of
# magnitude, while the quantile-based credible interval (which the forecast
# actually reports) is unaffected. Winsorizing yields a spread consistent with
# that interval. The default trim = 0.01 clips the outer 1% of draws per side:
# enough to remove the single-draw domination of a spurious near-non-stationary
# tail (e.g. one |phi| >= 1 draw in a mostly-stationary posterior), while leaving
# genuine non-stationary posterior mass visible (a species whose AR posterior is
# truly explosive still gets a large temporal variance). For a well-behaved
# column the winsorized variance is within ~5% of the raw variance.
.robust_col_var <- function(mat, trim = 0.01) {
  if (trim <= 0) return(apply(mat, 2, stats::var))
  apply(mat, 2, function(col) {
    qs <- stats::quantile(col, probs = c(trim, 1 - trim),
                          names = FALSE, na.rm = TRUE)
    stats::var(pmin(pmax(col, qs[1]), qs[2]))
  })
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
# All components are computed on the RESPONSE SCALE (g^{-1}(eta)) so that the
# law of total variance holds and the components are dimensionally consistent.
# total_var is DEFINED as the sum of its components, so the budget reconciles
# exactly and the percentage shares sum to 100%:
#   iid       : V_total = V_param + V_env + V_resid
#   has_autocor: V_total = V_param + V_env + V_resid + V_temporal
# V_env (perturbed predictors) and V_temporal (AR/MA/ARMA accumulation, the
# excess of the sampled posterior-predictive variance over V_param + V_resid)
# are therefore genuine sub-shares of a total that contains them, not
# quantities bolted on top of an env-free total.
#
# pp           : [n_draws  x n_obs] posterior predictive draws (response scale)
# mu_draws     : [n_draws  x n_obs] response-scale posterior means
#                (= posterior_linpred(transform=TRUE); equals lp for Gaussian)
# mu_perturbed : [n_perturb x n_obs] response-scale draws under env perturbation
# mu_draws_sub : [n_perturb x n_obs] first n_perturb rows of mu_draws
# family       : brms family object (used to dispatch residual variance)
# disp_draws   : named list of dispersion parameter draws (sigma, phi, nu, shape)
# has_autocor  : logical; if TRUE, add a temporal_var component capturing the
#                AR/MA/ARMA-induced predictive variance beyond the iid sum.
.decompose_from_arrays <- function(pp, mu_draws, mu_perturbed, mu_draws_sub,
                                    family, disp_draws,
                                    has_autocor = FALSE) {
  n_obs <- ncol(pp)
  n_p   <- nrow(mu_perturbed)

  # Parameter uncertainty: variance of the response-scale posterior mean.
  param_var <- apply(mu_draws, 2, stats::var)

  # Environmental variance: subtraction estimator on the response scale.
  # Both column variances must be computed the SAME way: stats::var() dispatches
  # on na.rm to a different C code path in cov.c ("na.or.complete" vs
  # "everything"), and those two paths are not guaranteed to be bit-identical on
  # every platform (they are on x86_64, but not on aarch64). An asymmetric call
  # would make the difference below nonzero in the last bits even when the two
  # matrices are the same sample.
  v_perturbed <- apply(mu_perturbed, 2, stats::var, na.rm = TRUE)
  v_param_sub <- apply(mu_draws_sub, 2, stats::var, na.rm = TRUE)
  v_env_raw   <- v_perturbed - v_param_sub

  # When there is no environmental noise the two matrices are the same sample
  # (all_zero_noise path), so v_env_raw is 0 in exact arithmetic. Floating-point
  # summation is not bit-reproducible across platforms/BLAS, so test against a
  # relative tolerance rather than == 0. The threshold (sqrt(machine epsilon),
  # the all.equal() convention) is many orders of magnitude below anything the
  # subtraction estimator could resolve: its own Monte Carlo error is
  # O(var * sqrt(2 / n_p)), so an env share this small is numerically
  # indistinguishable from zero and is reported as zero.
  zero_tol <- sqrt(.Machine$double.eps) * pmax(v_perturbed, v_param_sub, 0)
  is_zero  <- !is.na(v_env_raw) & abs(v_env_raw) <= zero_tol

  env_var           <- pmax(v_env_raw, 0)
  env_var[is_zero]  <- 0

  # Monte Carlo SE of env_var.
  # Uses the chi-squared approximation SE(Var) ~= Var * sqrt(2 / (n - 1)).
  # Exactly zero on the no-noise path (env_var = 0 by construction); still
  # reported when v_env_raw is negative, since there the SE is what tells the
  # user the estimate is Monte Carlo noise.
  f_se       <- sqrt(2 / max(n_p - 1, 1))
  v_env_mcse <- ifelse(
    is_zero,
    0,
    sqrt((v_perturbed * f_se)^2 + (v_param_sub * f_se)^2)
  )

  # Residual variance: family-specific expected within-draw variance.
  # For models with an autocor term this is the *innovation* variance
  # (E[sigma^2] for Gaussian AR(1)); the additional autocorrelation-induced
  # spread is captured in temporal_var below.
  residual_var <- .family_residual_var(family, mu_draws, disp_draws)

  # Sampled posterior-predictive variance on the UNPERTURBED predictors, taken
  # tail-robustly (winsorized). By the law of total variance this equals
  # param_var + residual_var for iid models (up to Monte Carlo error) and
  # additionally contains the AR/MA/ARMA accumulation for autocorrelation models
  # (brms forecasts iteratively). It does NOT contain env_var, which arises only
  # from perturbing the predictors. The winsorization matters only for the
  # autocorrelation branch: a near-non-stationary AR posterior makes the far-lead
  # forecast variance heavy-tailed, so a raw var() would be dominated by a few
  # |phi| >= 1 draws and disagree with the quantile credible interval; the robust
  # estimate stays consistent with the reported interval (see .robust_col_var).
  pp_var_raw <- apply(pp, 2, stats::var)
  pp_var     <- .robust_col_var(pp)

  if (isTRUE(has_autocor)) {
    # Temporal (autocorrelation-induced) variance: the excess of the (robust)
    # posterior-predictive variance over the iid-equivalent (param + residual)
    # sum. For Gaussian AR(p) this is the autocorrelated accumulation of
    # innovations beyond a single sigma^2; for ARMA / MA / cosy() / sar() /
    # car() it absorbs whatever residual dependence structure brms is
    # modelling — the part the analytic single-innovation residual_var cannot.
    # Clamped at 0 to absorb small Monte Carlo noise. env is excluded from
    # this gap because pp uses unperturbed predictors (so env is not in pp_var
    # and must not be double-subtracted); it re-enters total_var below.
    temporal_var <- pmax(0, pp_var - (param_var + residual_var))
    # Total variance now CONTAINS every component. env enters through the
    # perturbed predictors, temporal through the AR accumulation. Defining
    # total as the component sum makes the budget reconcile EXACTLY (law of
    # total variance), so the reported percentage shares sum to 100%.
    total_var <- param_var + env_var + residual_var + temporal_var
    # Heads-up when the AR forecast variance is heavy-tailed (near / above the
    # unit root): the raw predictive variance greatly exceeds its winsorized
    # value, so the point forecast becomes uninformative faster than the robust
    # interval suggests. The reported components use the robust estimate.
    # Compare the mean raw vs mean robust predictive variance (mean, not max,
    # so a single explosive draw in one row does not by itself trip the flag).
    tail_ratio <- mean(pp_var_raw, na.rm = TRUE) /
                  max(mean(pp_var, na.rm = TRUE), .Machine$double.eps)
    if (is.finite(tail_ratio) && tail_ratio > 3) {
      .et_warn("Autocorrelation forecast variance is heavy-tailed (raw / robust ",
               "predictive variance up to ", round(tail_ratio, 1), "x): the AR ",
               "coefficient's posterior has mass near or above the unit root, so ",
               "a few non-stationary draws inflate the raw variance. Reported ",
               "temporal_var / total_var use the tail-robust estimate, consistent ",
               "with the quantile credible interval; treat long-lead point ",
               "forecasts with caution and see et_diagnose().")
    }
  } else {
    temporal_var <- NULL
    # iid: total variance is the exact law-of-total-variance sum. Using the
    # analytic sum (rather than the noisier sampled pp variance) guarantees
    # param + env + residual == total exactly.
    total_var <- param_var + env_var + residual_var
  }

  out <- data.frame(
    obs_id       = seq_len(n_obs),
    total_var    = total_var,
    param_var    = param_var,
    env_var      = env_var,
    v_env_mcse   = v_env_mcse,
    residual_var = residual_var,
    stringsAsFactors = FALSE
  )
  if (!is.null(temporal_var)) out$temporal_var <- temporal_var

  out
}

# ******************************************************************************
# Internal: family-aware helpers
# ______________________________________________________________________________

# Returns TRUE for Gaussian family with identity link.
# For this case posterior_linpred == posterior_linpred(transform=TRUE) and
# the decomposition is identical to the pre-response-scale implementation.
.is_gaussian_identity <- function(family) {
  if (is.null(family)) return(TRUE)
  fname <- tolower(family$family %||% "")
  lname <- tolower(family$link   %||% "identity")
  fname %in% c("gaussian", "normal") && lname == "identity"
}

# Apply the inverse link function g^{-1} to a matrix of link-scale values.
# Returns a matrix of the same dimensions on the response scale.
# Falls back to the identity (no transform) when the link is unknown.
.apply_inv_link <- function(lp_mat, family) {
  if (is.null(family)) return(lp_mat)
  link_str <- family$link %||% "identity"
  if (link_str == "identity") return(lp_mat)

  linkinv <- tryCatch(
    family$linkinv %||% stats::make.link(link_str)$linkinv,
    error = function(e) NULL
  )
  if (is.null(linkinv)) {
    linkinv <- tryCatch(stats::make.link(link_str)$linkinv,
                        error = function(e) NULL)
  }
  if (is.null(linkinv)) {
    .et_warn("Could not determine inverse link for '", link_str,
             "'; response-scale decomposition may be on the link scale.")
    return(lp_mat)
  }

  matrix(linkinv(as.vector(lp_mat)), nrow = nrow(lp_mat), ncol = ncol(lp_mat))
}

# Compute the expected within-draw variance E_s[Var(y | mu^(s))] for the
# supported brms families. Returns a numeric vector of length n_obs.
#
# For Gaussian the residual variance is constant across observations (sigma^2);
# for all other supported families it is observation-specific because the
# variance function V(mu) depends on the fitted mean.
#
# Supported families and their variance functions V(mu):
#   gaussian / normal : sigma^2                        (scalar, constant)
#   student           : sigma^2 * nu / (nu - 2)       (scalar, nu > 2)
#   binomial / bernoulli : mu * (1 - mu)
#   poisson           : mu
#   negbinomial       : mu + mu^2 / shape
#   beta              : mu * (1 - mu) / (phi + 1)
#   gamma             : mu^2 / shape
#
# Unsupported families return NA with a warning.
.family_residual_var <- function(family, mu_draws, disp_draws) {
  n_obs <- ncol(mu_draws)
  n_s   <- nrow(mu_draws)
  fname <- tolower(family$family %||% "gaussian")

  # Helper: fetch first n_s rows of a named dispersion vector.
  get_d <- function(name) {
    v <- disp_draws[[name]]
    if (is.null(v) || all(is.na(v))) return(NULL)
    v[seq_len(min(n_s, length(v)))]
  }

  if (fname %in% c("gaussian", "normal")) {
    sigma <- get_d("sigma")
    v <- if (!is.null(sigma)) mean(sigma^2, na.rm = TRUE) else NA_real_
    return(rep(v, n_obs))
  }

  if (fname == "student") {
    sigma <- get_d("sigma")
    nu    <- get_d("nu")
    if (is.null(sigma) || is.null(nu)) return(rep(NA_real_, n_obs))
    valid <- !is.na(nu) & nu > 2
    v <- if (any(valid)) {
      mean(sigma[valid]^2 * nu[valid] / (nu[valid] - 2), na.rm = TRUE)
    } else NA_real_
    return(rep(v, n_obs))
  }

  if (fname %in% c("binomial", "bernoulli")) {
    # V(mu) = mu * (1 - mu), averaged over posterior draws.
    return(colMeans(mu_draws * (1 - mu_draws), na.rm = TRUE))
  }

  if (fname == "poisson") {
    # V(mu) = mu (mean == variance for Poisson).
    return(colMeans(mu_draws, na.rm = TRUE))
  }

  if (fname %in% c("negbinomial", "negbinomial2", "neg_binomial_2")) {
    # V(mu, shape) = mu + mu^2 / shape.
    shape <- get_d("shape")
    if (is.null(shape)) {
      .et_warn("negbinomial 'shape' draws not found; falling back to Poisson residual.")
      return(colMeans(mu_draws, na.rm = TRUE))
    }
    # shape is length n_s; broadcast over n_obs columns via sweep.
    return(colMeans(
      mu_draws + sweep(mu_draws^2, 1, shape, "/"),
      na.rm = TRUE
    ))
  }

  if (fname == "beta") {
    # V(mu, phi) = mu * (1 - mu) / (phi + 1) — brms uses precision phi.
    phi <- get_d("phi")
    if (is.null(phi)) {
      .et_warn("beta 'phi' draws not found; returning mu*(1-mu) approximation.")
      return(colMeans(mu_draws * (1 - mu_draws), na.rm = TRUE))
    }
    return(colMeans(
      sweep(mu_draws * (1 - mu_draws), 1, 1 + phi, "/"),
      na.rm = TRUE
    ))
  }

  if (fname %in% c("gamma", "gamma2")) {
    # V(mu, shape) = mu^2 / shape.
    shape <- get_d("shape")
    if (is.null(shape)) return(rep(NA_real_, n_obs))
    return(colMeans(
      sweep(mu_draws^2, 1, shape, "/"),
      na.rm = TRUE
    ))
  }

  .et_warn(
    "Family '", fname, "' is not supported for the response-scale residual ",
    "variance decomposition. residual_var will be NA. ",
    "Supported: gaussian, student, binomial, bernoulli, poisson, ",
    "negbinomial, beta, gamma."
  )
  rep(NA_real_, n_obs)
}

# Extract dispersion parameter draws from the full draws matrix.
# Always returns a named list with elements sigma, phi, shape, nu.
# Missing columns produce NULL entries.
.extract_disp_draws <- function(draws_mat, n_draws) {
  cols    <- colnames(draws_mat)
  n_use   <- min(n_draws, nrow(draws_mat))
  get_col <- function(name) {
    if (name %in% cols) draws_mat[seq_len(n_use), name] else NULL
  }
  list(
    sigma = get_col("sigma"),
    phi   = get_col("phi"),
    shape = get_col("shape"),
    nu    = get_col("nu")
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
  cat(sprintf("    Parameter    : %.4f\n", mean(decomp$param_var)))
  cat(sprintf("    Environmental: %.4f  (MC SE: %.4f)\n",
              mean(decomp$env_var), mean(decomp$v_env_mcse)))
  cat(sprintf("    Residual     : %.4f\n", mean(decomp$residual_var, na.rm = TRUE)))
  if ("temporal_var" %in% colnames(decomp)) {
    cat(sprintf("    Temporal     : %.4f\n", mean(decomp$temporal_var, na.rm = TRUE)))
  }
  cat(sprintf("    Total        : %.4f\n", mean(decomp$total_var)))
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
