# R/skill.R -- et_skill_score(): null-model-relative forecast skill.
#
# Shelf life (shelf_life.R) measures forecast PRECISION on the response scale.
# It is deliberately a complement to an ACCURACY-relative forecast limit: a
# precise forecast can still be wrong. Following the ecological-forecasting
# convention (Petchey et al. 2015; Wesselkamp et al. 2025; the NEON challenge),
# et_skill_score() scores the model against a sensible NULL model -- a random
# walk (persistence) or climatology -- with the continuous ranked probability
# score (CRPS), and reports the null-relative forecast limit: the lead time at
# which the model stops beating the null. Coupling the two (see the `skill`
# gate in shelf_life()) prevents a precise-but-biased forecast from passing.

#' Null-model-relative forecast skill (CRPS)
#'
#' Scores a forecast against a null model with the continuous ranked
#' probability score (CRPS) and returns the per-lead-time skill,
#' \eqn{1 - \mathrm{CRPS_{model}} / \mathrm{CRPS_{null}}} (positive when the
#' model beats the null).  This is the accuracy-relative counterpart to
#' \code{\link{shelf_life}}'s precision criterion; report both.  The
#' \strong{forecast limit} --- the first lead time beginning a sustained run of
#' non-skillful periods --- is attached as an attribute.
#'
#' @param predictions An \code{et_prediction} or \code{et_prediction_list}.
#' @param observed A \code{data.frame} of true responses positionally matched to
#'   \code{predictions$newdata} (with the grouping column for grouped input).
#' @param response_col Character.  Response column name; inferred from the model
#'   formula when \code{NULL}.
#' @param time_col Character.  Optional time column in \code{newdata} used for
#'   the \code{time} axis (lead time otherwise defaults to the row index).
#' @param null Character.  The null model: \code{"climatology"} (default; the
#'   unconditional response distribution) or \code{"random_walk"} (persistence:
#'   the last observed value carried forward with variance growing linearly with
#'   lead time).
#' @param climatology Optional numeric vector giving the climatological
#'   reference distribution.  Defaults to the training response values.
#' @param rw_sd,rw_start Optional random-walk per-step innovation SD and start
#'   value.  Default to \code{sd(diff(train_response))} and the last training
#'   response, respectively.
#' @param min_run Integer.  Number of \emph{consecutive} non-skillful periods
#'   required before the forecast limit is declared (default \code{2}).  Skill
#'   near zero crosses by chance, and a bare first-crossing rule
#'   (\code{min_run = 1}, the previous behaviour) promotes a lone dip to a
#'   horizon --- the same failure mode \code{\link{shelf_life}} guards against.
#'   Isolated dips are still reported via \code{n_below} / \code{first_below}.
#' @param ... Unused.
#' @return A \code{data.frame} with columns \code{obs_id}, \code{time},
#'   \code{crps_model}, \code{crps_null}, \code{skill}, and \code{skillful}
#'   (\code{skill > 0}); a leading \code{group} column for grouped input.  The
#'   attribute \code{"forecast_limit"} holds a list with the null-relative
#'   horizon, its \code{type} (\code{"observed"} or \code{"lower_bound"}), and
#'   the diagnostics \code{n_below}, \code{first_below} and \code{min_run};
#'   or is per-group for grouped input.
#' @seealso \code{\link{shelf_life}}, \code{\link{et_calibrate}}
#' @export
et_skill_score <- function(predictions, observed, response_col = NULL,
                           time_col = NULL,
                           null = c("climatology", "random_walk"),
                           climatology = NULL, rw_sd = NULL, rw_start = NULL,
                           min_run = 2L, ...) {
  UseMethod("et_skill_score")
}

#' @export
et_skill_score.et_prediction <- function(predictions, observed,
                                          response_col = NULL, time_col = NULL,
                                          null = c("climatology", "random_walk"),
                                          climatology = NULL, rw_sd = NULL,
                                          rw_start = NULL, min_run = 2L, ...) {
  null <- match.arg(null)
  response_col <- .resolve_response_col(response_col, predictions$model)
  if (nrow(observed) != nrow(predictions$newdata)) {
    stop("observed must have the same number of rows as newdata (",
         nrow(predictions$newdata), "), got ", nrow(observed), ".")
  }
  y_true <- observed[[response_col]]
  if (is.null(y_true)) {
    stop("Column '", response_col, "' not found in observed data.frame.")
  }
  train_resp <- .training_response(predictions$model, response_col)

  res <- .compute_skill(
    draws       = predictions$predictive_draws %||% predictions$posterior_predict,
    y_true      = y_true,
    newdata     = predictions$newdata,
    time_col    = time_col,
    null        = null,
    train_resp  = train_resp,
    climatology = climatology,
    rw_sd       = rw_sd,
    rw_start    = rw_start,
    min_run     = min_run
  )
  structure(res$df, forecast_limit = res$limit, null_model = null)
}

#' @export
et_skill_score.et_prediction_list <- function(predictions, observed,
                                               response_col = NULL,
                                               time_col = NULL,
                                               null = c("climatology", "random_walk"),
                                               climatology = NULL, rw_sd = NULL,
                                               rw_start = NULL, min_run = 2L, ...) {
  null <- match.arg(null)
  grouping <- predictions$grouping
  if (!grouping %in% colnames(observed)) {
    stop("Grouping column '", grouping, "' not found in observed.")
  }
  limit_by_group <- list()
  parts <- lapply(names(predictions$predictions), function(g) {
    pred <- predictions$predictions[[g]]
    if (is.null(pred)) return(NULL)
    obs_g <- observed[observed[[grouping]] == g, , drop = FALSE]
    if (nrow(obs_g) == 0) return(NULL)
    rc <- .resolve_response_col(response_col, pred$model)
    y_true <- obs_g[[rc]]
    if (is.null(y_true) || length(y_true) != nrow(pred$newdata)) return(NULL)
    res <- .compute_skill(
      draws       = pred$predictive_draws %||% pred$posterior_predict,
      y_true      = y_true, newdata = pred$newdata, time_col = time_col,
      null = null, train_resp = .training_response(pred$model, rc),
      climatology = climatology, rw_sd = rw_sd, rw_start = rw_start,
      min_run = min_run)
    limit_by_group[[g]] <<- res$limit
    cbind(data.frame(group = g, stringsAsFactors = FALSE), res$df)
  })
  out <- do.call(rbind, Filter(Negate(is.null), parts))
  structure(out, forecast_limit_by_group = limit_by_group, null_model = null)
}

# --- internals ---------------------------------------------------------------

.training_response <- function(model, response_col) {
  if (is.null(model) || is.null(model$data) ||
      !response_col %in% colnames(model$data)) {
    return(NULL)
  }
  as.numeric(model$data[[response_col]])
}

.compute_skill <- function(draws, y_true, newdata, time_col, null,
                           train_resp, climatology, rw_sd, rw_start,
                           min_run = 2L) {
  n <- length(y_true)
  time_vals <- if (!is.null(time_col) && time_col %in% colnames(newdata)) {
    newdata[[time_col]]
  } else {
    seq_len(n)
  }

  crps_model <- vapply(seq_len(n), function(i) {
    yi <- y_true[i]
    if (is.na(yi)) return(NA_real_)
    col <- draws[, i]; col <- col[!is.na(col)]
    if (length(col) == 0L) return(NA_real_) else .crps_sample(col, yi)
  }, numeric(1))

  if (null == "climatology") {
    ref <- climatology %||% train_resp
    if (is.null(ref) || length(ref) < 2L) {
      stop("climatology null needs a reference distribution; supply ",
           "`climatology=` or ensure the model carries training data.")
    }
    crps_null <- vapply(y_true, function(yi)
      if (is.na(yi)) NA_real_ else .crps_sample(ref, yi), numeric(1))
  } else {                                   # random_walk (persistence)
    if (is.null(rw_start)) {
      if (is.null(train_resp)) stop("random_walk null needs training data or rw_start.")
      rw_start <- train_resp[length(train_resp)]
    }
    if (is.null(rw_sd)) {
      if (is.null(train_resp) || length(train_resp) < 3L)
        stop("random_walk null needs training data or rw_sd.")
      rw_sd <- stats::sd(diff(train_resp))
    }
    lead <- seq_len(n)                       # one RW step per forecast step
    sigma_k <- rw_sd * sqrt(lead)
    crps_null <- vapply(seq_len(n), function(i)
      if (is.na(y_true[i])) NA_real_
      else .crps_norm(y_true[i], rw_start, sigma_k[i]), numeric(1))
  }

  skill <- 1 - crps_model / crps_null
  df <- data.frame(
    obs_id = seq_len(n), time = time_vals,
    crps_model = crps_model, crps_null = crps_null,
    skill = skill, skillful = skill > 0,
    stringsAsFactors = FALSE
  )

  # Forecast limit: the first lead time beginning a SUSTAINED run of
  # non-skillful periods. A bare first-crossing rule has the same failure mode
  # here that it has in shelf_life(): when skill hovers near zero without
  # trending, one lead dips below by chance and is promoted to "the horizon".
  # Kyoto against climatology is a live example -- mean skill is comfortably
  # positive across the window, yet a single negative lead would set the limit
  # at the seventh forecast year. Isolated dips are still reported, as
  # diagnostics, so they are visible rather than silently load-bearing.
  not_ok      <- !df$skillful & !is.na(df$skill)
  first_below <- if (any(not_ok)) df$time[which(not_ok)[1]] else NA
  idx         <- .first_sustained_index(not_ok, min_run)
  limit <- if (!is.na(idx)) {
    list(value = df$time[idx], type = "observed",
         n_below = sum(not_ok), first_below = first_below, min_run = min_run,
         description = paste0("Model stops beating the ", null,
                              " null at time ", df$time[idx],
                              " (sustained over ", min_run, " periods)."))
  } else {
    list(value = NA_real_, type = "lower_bound",
         n_below = sum(not_ok), first_below = first_below, min_run = min_run,
         description = paste0("Model beats the ", null,
                              " null throughout (limit > ", max(df$time),
                              "); ", sum(not_ok), " isolated non-skillful ",
                              "period(s) did not persist for ", min_run, "."))
  }
  list(df = df, limit = limit)
}

# CRPS of an empirical sample x against observation y (O(S log S) estimator).
.crps_sample <- function(x, y) {
  S  <- length(x)
  xs <- sort(x)
  mean(abs(x - y)) - sum((2 * seq_len(S) - S - 1) * xs) / S^2
}

# CRPS of a Normal(mu, sigma) predictive against observation y (closed form).
.crps_norm <- function(y, mu, sigma) {
  if (sigma <= 0) return(abs(y - mu))
  z <- (y - mu) / sigma
  sigma * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
}
