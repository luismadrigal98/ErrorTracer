# R/shelf_life.R -- shelf_life(): forecast horizon / shelf life analysis

#' Compute the forecast shelf life
#'
#' Quantifies \emph{when} a forecast becomes uninformative by comparing the
#' width of credible intervals to a response scale.  A forecast is
#' uninformative when its CI width exceeds \code{threshold * response_scale}.
#'
#' Shelf life is a \strong{precision} criterion on the response scale --- a
#' deliberate \emph{complement} to accuracy-relative forecast limits (Petchey
#' et al. 2015; Wesselkamp et al. 2025), not a replacement.  It answers "when do
#' the error bars get too wide to act on?", whereas a forecast limit answers
#' "when does the model stop beating a null?".  Because a forecast can be
#' precise yet biased, pass a \code{skill} table from \code{\link{et_skill_score}}
#' to \strong{gate} the shelf life on skill: a period then counts as
#' informative only if it is both precise \emph{and} beats the null model.  This
#' couples the two views and prevents a precise-but-biased forecast from
#' silently passing.
#'
#' The function operates in three modes depending on the available data and
#' whether the uninformative threshold is crossed within the forecast window:
#'
#' \describe{
#'   \item{\strong{Observed}}{The threshold is crossed within the
#'     forecast/validation window \emph{and the crossing persists}.  The
#'     shelf life is the first time point beginning a run of \code{min_run}
#'     consecutive periods with \code{ratio >= threshold}.  Requiring
#'     persistence matters: when the ratio hovers near \code{threshold}
#'     without trending, a single noisy period can exceed it by chance, and
#'     the first-exceedance rule used before ErrorTracer 1.3.1 would report
#'     that excursion as a horizon.  Isolated exceedances are now reported in
#'     the horizon diagnostics (\code{n_exceedances},
#'     \code{first_exceedance}) instead of being promoted to a crossing.}
#'   \item{\strong{Projected}}{No sustained crossing occurs in-window but the
#'     CI/range ratio is trending upward \emph{detectably} --- the fitted
#'     slope must exceed \code{min_slope_for_projection} \strong{and} be
#'     significantly positive at \code{projection_alpha}.  A linear trend is
#'     extrapolated to estimate when the threshold would be reached.  The
#'     projected crossing time \eqn{t^* = (\tau - a) / b} (where \eqn{\tau}
#'     is the threshold, \eqn{a} the fitted intercept, \eqn{b} the fitted
#'     slope) is reported together with a standard error \code{se_t_star}
#'     derived via the delta method.}
#'   \item{\strong{Lower bound}}{No sustained crossing and no trend
#'     distinguishable from flat (or a projection beyond the extrapolation
#'     cap).  The shelf life is reported as a lower bound:
#'     \code{> last observed time}.  This is the honest outcome for a
#'     forecast whose precision simply does not degrade in-window.}
#' }
#'
#' The intended workflow is:
#' \enumerate{
#'   \item Fit the model (\code{\link{et_fit}}).
#'   \item Predict on held-out data or a future time window
#'     (\code{\link{et_predict}}).
#'   \item Call \code{shelf_life()} on those predictions.
#'   \item If the threshold is not crossed in the held-out window, the
#'     \emph{projected} mode extrapolates the horizon automatically.
#' }
#'
#' @param predictions An \code{et_prediction} or \code{et_prediction_list}.
#' @param response_scale Numeric vector of length 2 (\code{c(min, max)})
#'   giving the response scale used as the denominator in the
#'   CI-width / range ratio.  Supply a biologically motivated interval when
#'   one exists.  When \code{NULL} (default), the observed range of the
#'   \strong{training response} is used as a documented, reproducible default
#'   (with an informational message) instead of an arbitrary number.
#'   The effective range is \code{diff(response_scale)}.
#' @param ci_level Numeric.  The credible interval level to use (default
#'   0.90).  Must be present in the \code{et_prediction} object.
#' @param threshold Numeric.  CI width / response scale above which the
#'   forecast is uninformative (default 1.0).
#' @param time_col Character.  Optional column in
#'   \code{predictions$newdata} to use as the time axis.  If \code{NULL},
#'   observation index is used.
#' @param skill Optional skill gate.  A \code{data.frame} returned by
#'   \code{\link{et_skill_score}} (positionally matched to the forecast rows;
#'   its \code{skillful} column is used) or a bare logical vector.  When
#'   supplied, a forecast period is \code{informative} only if it is both
#'   precise \emph{and} skillful, so a precise-but-biased forecast does not
#'   pass.  When \code{NULL} (default) the precision criterion is used alone.
#' @param min_slope_for_projection Numeric.  Minimum linear slope (of
#'   ratio vs.\ time) required to attempt extrapolation when all periods
#'   are informative.  Below this value the shelf life is reported as a
#'   lower bound only.  Default \code{1e-4}.
#' @param max_extrapolation_factor Numeric.  Cap on how far the linear
#'   projection may reach beyond the observed window.  If the projected
#'   crossing time exceeds
#'   \code{max(time) + max_extrapolation_factor * (max(time) - min(time))},
#'   the result is reported as a lower bound instead of a projection.
#'   Set to \code{Inf} to disable the cap.  Default \code{10}.
#' @param min_run Integer.  Number of \emph{consecutive} uninformative
#'   periods required before the threshold counts as crossed (default
#'   \code{2}).  Guards against a single noisy period manufacturing a
#'   horizon when the ratio hovers near \code{threshold}.  Set
#'   \code{min_run = 1} to recover the pre-1.3.1 first-exceedance
#'   behaviour; raise it for noisier series.
#' @param projection_alpha Numeric.  Significance level for the one-sided
#'   test that the fitted ratio-vs-time slope is positive (default
#'   \code{0.05}).  A slope can clear \code{min_slope_for_projection} while
#'   being statistically indistinguishable from flat; extrapolating it
#'   yields a confident-looking crossing time whose interval may span
#'   centuries and include the past.  When the slope fails this test the
#'   result is reported as a lower bound instead of a projection.  Set to
#'   \code{1} to disable the test and recover the pre-1.3.1 behaviour.
#' @param ... Unused.
#' @param plausible_range Deprecated.  Use \code{response_scale} instead.
#'
#' @return An \code{et_shelf_life} object (a \code{data.frame}) with columns:
#'   \describe{
#'     \item{obs_id}{Observation index.}
#'     \item{time}{Time axis value.}
#'     \item{ci_width}{Width of the credible interval at \code{ci_level}.}
#'     \item{plausible_range}{Effective response scale (scalar diff).}
#'     \item{ratio}{CI width / response scale.}
#'     \item{informative}{Logical; \code{TRUE} when ratio < threshold.}
#'   }
#'   For grouped predictions a \code{group} column is prepended.
#'   The object carries three attributes that drive \code{print()}:
#'   \describe{
#'     \item{\code{horizon}}{Named list with elements \code{value},
#'       \code{type} (\code{"observed"}, \code{"projected"}, or
#'       \code{"lower_bound"}), \code{last_informative},
#'       \code{description}, (for projected) \code{se_t_star}, and the
#'       exceedance diagnostics \code{n_exceedances},
#'       \code{frac_exceedance}, \code{first_exceedance} and
#'       \code{min_run}.  Comparing \code{first_exceedance} with
#'       \code{value} shows whether the reported horizon rests on a
#'       persistent crossing or whether isolated excursions preceded it.}
#'     \item{\code{horizon_by_group}}{For grouped objects: named list of
#'       per-group horizon lists.}
#'     \item{\code{threshold}}{The threshold value used.}
#'   }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' df  <- data.frame(y = rnorm(20), year = 2001:2020, x1 = rnorm(20))
#' fit <- et_fit(y ~ x1, data = df,
#'               chains = 1, iter = 500, warmup = 250,
#'               cores = 1, refresh = 0)
#' new_df <- data.frame(x1 = rnorm(10), year = 2021:2030)
#' pred   <- et_predict(fit, newdata = new_df,
#'                      n_draws = 200, n_perturb = 50)
#' sl <- shelf_life(pred,
#'                  response_scale            = range(df$y),
#'                  ci_level                  = 0.90,
#'                  threshold                 = 1.0,
#'                  time_col                  = "year",
#'                  max_extrapolation_factor  = 10)
#' print(sl)
#' }
#' @export
shelf_life <- function(predictions,
                       response_scale = NULL,
                       ci_level = 0.90,
                       threshold = 1.0,
                       time_col = NULL,
                       skill = NULL,
                       min_slope_for_projection = 1e-4,
                       max_extrapolation_factor = 10,
                       min_run = 2L,
                       projection_alpha = 0.05,
                       ...,
                       plausible_range = NULL) {
  # Deprecated alias: plausible_range -> response_scale. Re-call (rather than
  # reassign) because UseMethod() dispatches with the originally-supplied
  # arguments, not values modified in this frame.
  if (!is.null(plausible_range)) {
    if (!is.null(response_scale)) {
      stop("Cannot specify both 'response_scale' and 'plausible_range'.")
    }
    warning(
      "'plausible_range' is deprecated; use 'response_scale' instead.",
      call. = FALSE
    )
    return(shelf_life(
      predictions              = predictions,
      response_scale           = plausible_range,
      ci_level                 = ci_level,
      threshold                = threshold,
      time_col                 = time_col,
      skill                    = skill,
      min_slope_for_projection = min_slope_for_projection,
      max_extrapolation_factor = max_extrapolation_factor,
      min_run                  = min_run,
      projection_alpha         = projection_alpha,
      ...
    ))
  }
  UseMethod("shelf_life")
}

#' @export
shelf_life.et_prediction <- function(predictions,
                                     response_scale = NULL,
                                     ci_level = 0.90,
                                     threshold = 1.0,
                                     time_col = NULL,
                                     skill = NULL,
                                     min_slope_for_projection = 1e-4,
                                     max_extrapolation_factor = 10,
                                     min_run = 2L,
                                     projection_alpha = 0.05,
                                     ...,
                                     plausible_range = NULL) {
  skill_ok <- if (!is.null(skill)) .skill_ok_vector(skill) else NULL
  .compute_shelf_life_single(
    predictions = predictions,
    response_scale = response_scale,
    ci_level = ci_level,
    threshold = threshold,
    time_col = time_col,
    skill_ok = skill_ok,
    min_slope_for_projection = min_slope_for_projection,
    max_extrapolation_factor = max_extrapolation_factor,
    min_run = min_run,
    projection_alpha = projection_alpha
  )
}

#' @export
shelf_life.et_prediction_list <- function(predictions,
                                          response_scale = NULL,
                                          ci_level = 0.90,
                                          threshold = 1.0,
                                          time_col = NULL,
                                          skill = NULL,
                                          min_slope_for_projection = 1e-4,
                                          max_extrapolation_factor = 10,
                                          min_run = 2L,
                                          projection_alpha = 0.05,
                                          ...,
                                          plausible_range = NULL) {
  horizon_by_group <- list()

  parts <- lapply(names(predictions$predictions), function(g) {
    pred <- predictions$predictions[[g]]
    if (is.null(pred)) return(NULL)

    # Per-group skill gate: subset the skill table to this group if grouped.
    skill_ok_g <- if (!is.null(skill)) {
      if ("group" %in% colnames(skill)) {
        .skill_ok_vector(skill[skill$group == g, , drop = FALSE])
      } else {
        .skill_ok_vector(skill)
      }
    } else NULL

    sl <- .compute_shelf_life_single(
      pred, response_scale, ci_level, threshold,
      time_col, skill_ok = skill_ok_g,
      min_slope_for_projection = min_slope_for_projection,
      max_extrapolation_factor = max_extrapolation_factor,
      min_run = min_run,
      projection_alpha = projection_alpha
    )
    horizon_by_group[[g]] <<- attr(sl, "horizon")
    cbind(data.frame(group = g, stringsAsFactors = FALSE), sl)
  })

  result <- do.call(rbind, Filter(Negate(is.null), parts))
  result <- structure(result, class = c("et_shelf_life", "data.frame"))
  attr(result, "horizon_by_group") <- horizon_by_group
  attr(result, "threshold") <- threshold
  result
}

#' @export
shelf_life.default <- function(predictions, ...) {
  stop("shelf_life() expects an et_prediction or et_prediction_list object.")
}

# ******************************************************************************
# Internal: compute shelf life for a single et_prediction
# ______________________________________________________________________________

.compute_shelf_life_single <- function(predictions, response_scale,
                                        ci_level, threshold, time_col,
                                        skill_ok = NULL,
                                        min_slope_for_projection,
                                        max_extrapolation_factor = 10,
                                        min_run = 2L,
                                        projection_alpha = 0.05) {
  ci_df <- predictions$credible_intervals
  avail_levels <- unique(ci_df$ci_level)

  if (!ci_level %in% avail_levels) {
    stop(
      "ci_level ", ci_level, " not found in predictions. ",
      "Available: ", paste(avail_levels, collapse = ", "),
      ". Re-run et_predict() with the desired level in ci_levels."
    )
  }

  ci_sub <- ci_df[ci_df$ci_level == ci_level, ]

  # Principled default: when response_scale is unspecified, use the observed
  # range of the training response as a conservative, documented default rather
  # than an arbitrary user number (addresses the "where did the threshold come
  # from?" critique). Users should still supply a domain-meaningful range when
  # one exists.
  if (is.null(response_scale)) {
    response_scale <- .default_response_scale(predictions)
    .et_info("shelf_life(): response_scale not supplied; using the training ",
             "response range c(", round(response_scale[1], 4), ", ",
             round(response_scale[2], 4), "). Supply a domain-meaningful ",
             "range for a principled threshold.")
  }
  if (length(response_scale) != 2) {
    stop("response_scale must be a numeric vector of length 2: c(min, max).")
  }
  pr <- abs(diff(response_scale))
  if (pr < 1e-12) stop("response_scale min and max are equal.")

  time_vals <- if (!is.null(time_col) &&
                   time_col %in% colnames(predictions$newdata)) {
    predictions$newdata[[time_col]]
  } else {
    ci_sub$row_id
  }

  # Precision criterion: CI width small relative to the response scale.
  precise <- ci_sub$width < threshold * pr
  # Skill gate: a precise forecast is only *informative* if it is also skillful
  # (beats the null / is not precise-but-biased). Without a skill table the gate
  # is a no-op (precision alone).
  if (!is.null(skill_ok)) {
    if (length(skill_ok) != length(precise)) {
      .et_warn("skill gate length (", length(skill_ok), ") does not match the ",
               "number of forecast rows (", length(precise),
               "); ignoring the skill gate.")
      informative <- precise
    } else {
      informative <- precise & (skill_ok %in% TRUE)
    }
  } else {
    informative <- precise
  }

  result_df <- data.frame(
    obs_id = ci_sub$row_id,
    time = time_vals,
    ci_width = ci_sub$width,
    plausible_range = pr,
    ratio = ci_sub$width / pr,
    informative = informative,
    stringsAsFactors = FALSE
  )

  # Compute shelf life horizon
  horizon <- .compute_horizon(result_df, threshold, min_slope_for_projection,
                              max_extrapolation_factor, min_run = min_run,
                              projection_alpha = projection_alpha)

  result <- structure(result_df, class = c("et_shelf_life", "data.frame"))
  attr(result, "horizon") <- horizon
  attr(result, "threshold") <- threshold
  result
}

# Extract the logical skill gate from an et_skill_score() table (or a bare
# logical vector).
.skill_ok_vector <- function(skill) {
  if (is.data.frame(skill) && "skillful" %in% colnames(skill)) {
    as.logical(skill$skillful)
  } else if (is.logical(skill)) {
    skill
  } else {
    stop("`skill` must be a data.frame from et_skill_score() (with a ",
         "'skillful' column) or a logical vector.")
  }
}

# Principled default response scale: the observed range of the training
# response. Documented and reproducible, unlike an arbitrary user number.
.default_response_scale <- function(predictions) {
  model <- predictions$model
  rc <- tryCatch(.resolve_response_col(NULL, model), error = function(e) NULL)
  if (!is.null(rc) && !is.null(model$data) && rc %in% colnames(model$data)) {
    rng <- range(model$data[[rc]], na.rm = TRUE)
    if (is.finite(diff(rng)) && diff(rng) > 1e-12) return(rng)
  }
  stop("response_scale not supplied and could not be derived from the model's ",
       "training data; please supply response_scale = c(min, max).")
}

# ******************************************************************************
# Internal: derive horizon from a single-group shelf-life table
# ______________________________________________________________________________


# First index at which `min_run` consecutive periods are all uninformative.
# Returns NA_integer_ when no such run exists. min_run = 1 reproduces the
# historical first-exceedance behaviour.
.first_sustained_index <- function(uninformative, min_run) {
  n <- length(uninformative)
  if (n == 0L) return(NA_integer_)
  if (min_run <= 1L) {
    return(if (any(uninformative)) which(uninformative)[1] else NA_integer_)
  }
  if (n < min_run) return(NA_integer_)
  for (i in seq_len(n - min_run + 1L)) {
    if (all(uninformative[i:(i + min_run - 1L)])) return(i)
  }
  NA_integer_
}

.compute_horizon <- function(sl_df, threshold, min_slope,
                             max_extrapolation_factor = 10,
                             min_run = 2L,
                             projection_alpha = 0.05) {
  times <- sl_df$time
  ratios <- sl_df$ratio
  informative <- sl_df$informative
  uninformative <- !informative

  # Diagnostics reported on every horizon, so an isolated excursion is visible
  # rather than silently promoted to a horizon (which the first-exceedance rule
  # used to do: a single noisy year manufactured a shelf life).
  n_exceed     <- sum(uninformative)
  first_exceed <- if (n_exceed > 0L) min(times[uninformative]) else NA_real_
  diag_fields  <- list(
    n_exceedances    = n_exceed,
    frac_exceedance  = n_exceed / length(uninformative),
    first_exceedance = first_exceed,
    min_run          = min_run
  )

  # --- Mode 1: threshold crossed and the crossing is SUSTAINED ---
  idx <- .first_sustained_index(uninformative, min_run)
  if (!is.na(idx)) {
    first_cross <- times[idx]
    last_ok     <- if (any(informative & seq_along(times) < idx)) {
      max(times[informative & seq_along(times) < idx])
    } else NA_real_
    return(c(list(
      value = first_cross,
      type = "observed",
      last_informative = last_ok,
      description = paste0(
        "Threshold (ratio >= ", threshold, ") sustained (>= ", min_run,
        " consecutive periods) from ", first_cross, "."
      )
    ), diag_fields))
  }

  # No sustained crossing. Any exceedances present are isolated excursions, so
  # fall through to the trend modes rather than reporting a spurious horizon.
  isolated_note <- if (n_exceed > 0L) {
    paste0(" ", n_exceed, " isolated exceedance(s) (first at ", first_exceed,
           ") did not persist for ", min_run,
           " consecutive periods and are not treated as a crossing.")
  } else ""

  # --- Mode 2: no sustained crossing -- try linear extrapolation ---
  last_ok <- max(times)

  if (is.numeric(times) && length(times) >= 3) {
    lm_fit <- suppressWarnings(lm(ratios ~ times))
    slope <- coef(lm_fit)[["times"]]
    b0 <- coef(lm_fit)[["(Intercept)"]]

    # The magnitude gate alone is not enough: a slope of 1e-3 clears it while
    # being statistically indistinguishable from flat, and extrapolating such a
    # slope yields a confident-looking t* whose interval can span centuries and
    # include the past. Require the slope to be significantly positive as well.
    slope_p <- tryCatch({
      cf <- suppressWarnings(summary(lm_fit))$coefficients
      stats::pt(cf["times", "t value"], df = lm_fit$df.residual,
                lower.tail = FALSE)          # one-sided: slope > 0
    }, error = function(e) NA_real_)
    slope_flat <- is.na(slope_p) || slope_p > projection_alpha

    if (!is.na(slope) && slope > min_slope && slope_flat) {
      return(c(list(
        value = NA_real_,
        type = "lower_bound",
        last_informative = last_ok,
        slope = slope,
        slope_p = slope_p,
        description = paste0(
          "No sustained exceedance across ", nrow(sl_df),
          " forecast periods. The CI/range ratio has no trend distinguishable ",
          "from flat (slope = ", round(slope, 5), ", one-sided p = ",
          signif(slope_p, 2), " > ", projection_alpha,
          "), so no crossing is projected. Shelf life > ", last_ok, ".",
          isolated_note
        )
      ), diag_fields))
    }

    if (!is.na(slope) && slope > min_slope) {
      proj <- (threshold - b0) / slope

      # SE of t* via delta method: t* = (tau - a) / b
      # d(t*)/d(a) = -1/b,  d(t*)/d(b) = -(tau - a)/b^2 = -t*/b
      vc   <- stats::vcov(lm_fit)
      va   <- vc["(Intercept)", "(Intercept)"]
      vb   <- vc["times", "times"]
      cab  <- vc["(Intercept)", "times"]
      se_t_star <- sqrt(max(0, (va + proj^2 * vb + 2 * proj * cab) / slope^2))

      # Cap extrapolation: if projection exceeds last_ok + factor * window,
      # the trend is too flat to give a meaningful estimate -- report lower bound.
      window   <- max(times) - min(times)
      cap_time <- max(times) + max_extrapolation_factor * window

      if (is.finite(max_extrapolation_factor) && proj > cap_time) {
        return(c(list(
          value = NA_real_,
          type = "lower_bound",
          last_informative = last_ok,
          slope = slope,
          slope_p = slope_p,
          description = paste0(
            "No sustained exceedance across ", nrow(sl_df),
            " forecast periods. Linear trend (slope = ", round(slope, 5),
            " per time unit) projects threshold crossing at ~", round(proj, 1),
            ", but this exceeds the extrapolation cap (",
            round(cap_time, 1), "). Shelf life > ", last_ok, ".",
            isolated_note
          )
        ), diag_fields))
      }

      return(c(list(
        value = proj,
        type = "projected",
        last_informative = last_ok,
        se_t_star = se_t_star,
        slope = slope,
        slope_p = slope_p,
        description = paste0(
          "No sustained exceedance across ", nrow(sl_df),
          " forecast periods. Linear trend (slope = ", round(slope, 5),
          " per time unit) projects threshold crossing at ~", round(proj, 1),
          " (SE = ", round(se_t_star, 2), ").", isolated_note
        )
      ), diag_fields))
    }
  }

  # --- Mode 3: no upward trend -- lower bound only ---
  c(list(
    value = NA_real_,
    type = "lower_bound",
    last_informative = last_ok,
    description = paste0(
      "No sustained exceedance across ", nrow(sl_df),
      " forecast periods and no upward trend in CI/range ratio. ",
      "Shelf life > ", last_ok, ".", isolated_note
    )
  ), diag_fields)
}

# ******************************************************************************
# Internal: format a horizon for printing
# ______________________________________________________________________________

.format_horizon <- function(hor) {
  if (is.null(hor)) return(NULL)
  switch(hor$type,
    observed  = paste0("~", hor$value,
                         " (observed -- threshold first exceeded)"),
    projected = paste0("~", round(hor$value, 1),
                         " (projected -- extrapolated beyond forecast window)"),
    lower_bound = paste0("> ", hor$last_informative,
                         " (lower bound -- all periods informative, no trend)")
  )
}

# ******************************************************************************
# S3 methods for et_shelf_life
# ______________________________________________________________________________

#' @export
print.et_shelf_life <- function(x, ...) {
  thr <- attr(x, "threshold")
  cat("ErrorTracer shelf life analysis\n")
  cat("  Observations     :", nrow(x), "\n")
  cat("  Plausible range  :", unique(x$plausible_range), "\n")
  cat("  Threshold        :", if (is.null(thr)) 1.0 else thr, "\n")

  if ("group" %in% colnames(x)) {
    groups  <- unique(x$group)
    hor_grp <- attr(x, "horizon_by_group")

    for (g in groups) {
      sub <- x[x$group == g, ]
      pct <- round(100 * mean(sub$informative), 1)
      hor <- if (!is.null(hor_grp)) hor_grp[[g]] else NULL
      hor_str <- if (!is.null(hor)) .format_horizon(hor) else ""
      cat(sprintf(
        "  [%s]  %.1f%% informative  mean ratio = %.3f  shelf life %s\n",
        g, pct, mean(sub$ratio), hor_str
      ))
    }
  } else {
    cat("  Informative      :", sum(x$informative), "/", nrow(x), "\n")
    cat("  Mean CI/range    :", round(mean(x$ratio), 3), "\n")
    cat("  Max CI/range     :", round(max(x$ratio), 3), "\n")

    hor <- attr(x, "horizon")
    if (!is.null(hor)) {
      cat("  Shelf life       :", .format_horizon(hor), "\n")
      cat("  Detail           :", hor$description, "\n")
    }
  }
  invisible(x)
}

#' @export
summary.et_shelf_life <- function(object, ...) {
  print(object)
  cat("\n--- Full table ---\n")
  print(as.data.frame(object))
  invisible(object)
}
