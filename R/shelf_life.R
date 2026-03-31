# R/shelf_life.R — shelf_life(): forecast horizon / shelf life analysis

#' Compute the forecast shelf life
#'
#' Quantifies \emph{when} forecasts become uninformative by comparing the
#' width of credible intervals to a plausible range of the response variable.
#' A forecast is considered uninformative when the CI width exceeds
#' \code{threshold * plausible_range}.
#'
#' This is the most novel concept in ErrorTracer and does not have a
#' direct equivalent in any existing R package.
#'
#' @param predictions An \code{et_prediction} or \code{et_prediction_list}.
#' @param plausible_range Numeric vector of length 2 (\code{c(min, max)})
#'   giving the biologically / ecologically plausible range of the response
#'   (e.g.\ \code{c(-1, 1)} for allele frequency change, or derived from
#'   training data).  The effective range is \code{diff(plausible_range)}.
#' @param ci_level Numeric.  The credible interval level to use for the
#'   shelf life calculation (default 0.90).  Must be present in the
#'   \code{et_prediction} object.
#' @param threshold Numeric.  The ratio (CI width / plausible range) above
#'   which the forecast is declared uninformative (default 1.0).
#'   A threshold of 1.0 means "CI spans the entire plausible range."
#' @param time_col Character.  Optional name of a column in
#'   \code{predictions$newdata} to use as the time axis label (e.g.\
#'   \code{"year"} or \code{"Period"}).  If \code{NULL}, the observation
#'   index is used.
#' @param ... Unused.
#' @return An \code{et_shelf_life} object (a \code{data.frame}) with columns:
#'   \describe{
#'     \item{obs_id}{Observation index (integer).}
#'     \item{time}{Time axis value (from \code{time_col} or obs_id).}
#'     \item{ci_width}{Width of the credible interval at \code{ci_level}.}
#'     \item{plausible_range}{The effective plausible range (scalar).}
#'     \item{ratio}{CI width / plausible range.}
#'     \item{informative}{Logical; \code{TRUE} when ratio < threshold.}
#'   }
#'   For grouped predictions, a \code{group} column is prepended.
#' @examples
#' \dontrun{
#' sl <- shelf_life(pred,
#'                  plausible_range = c(-1, 1),
#'                  ci_level = 0.95,
#'                  threshold = 1.0,
#'                  time_col = "year")
#' print(sl)
#' }
#' @export
shelf_life <- function(predictions, plausible_range, ci_level = 0.90,
                        threshold = 1.0, time_col = NULL, ...) {
  UseMethod("shelf_life")
}

#' @export
shelf_life.et_prediction <- function(predictions, plausible_range,
                                       ci_level = 0.90, threshold = 1.0,
                                       time_col = NULL, ...) {
  .compute_shelf_life_single(
    predictions    = predictions,
    plausible_range = plausible_range,
    ci_level       = ci_level,
    threshold      = threshold,
    time_col       = time_col
  )
}

#' @export
shelf_life.et_prediction_list <- function(predictions, plausible_range,
                                           ci_level = 0.90, threshold = 1.0,
                                           time_col = NULL, ...) {
  parts <- lapply(names(predictions$predictions), function(g) {
    pred <- predictions$predictions[[g]]
    if (is.null(pred)) return(NULL)
    sl <- .compute_shelf_life_single(pred, plausible_range, ci_level,
                                     threshold, time_col)
    cbind(data.frame(group = g, stringsAsFactors = FALSE), sl)
  })
  result <- do.call(rbind, Filter(Negate(is.null), parts))
  structure(result, class = c("et_shelf_life", "data.frame"))
}

#' @export
shelf_life.default <- function(predictions, ...) {
  stop("shelf_life() expects an et_prediction or et_prediction_list object.")
}

# ============================================================
# Internal
# ============================================================

.compute_shelf_life_single <- function(predictions, plausible_range,
                                        ci_level, threshold, time_col) {

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

  if (length(plausible_range) != 2) {
    stop("plausible_range must be a numeric vector of length 2: c(min, max).")
  }
  pr <- abs(diff(plausible_range))
  if (pr < 1e-12) stop("plausible_range min and max are equal.")

  time_vals <- if (!is.null(time_col) && time_col %in% colnames(predictions$newdata)) {
    predictions$newdata[[time_col]]
  } else {
    ci_sub$row_id
  }

  result <- data.frame(
    obs_id          = ci_sub$row_id,
    time            = time_vals,
    ci_width        = ci_sub$width,
    plausible_range = pr,
    ratio           = ci_sub$width / pr,
    informative     = ci_sub$width < threshold * pr,
    stringsAsFactors = FALSE
  )

  structure(result, class = c("et_shelf_life", "data.frame"))
}

# ============================================================
# S3 methods for et_shelf_life
# ============================================================

#' @export
print.et_shelf_life <- function(x, ...) {
  cat("ErrorTracer shelf life analysis\n")
  cat("  Observations     :", nrow(x), "\n")
  cat("  Plausible range  :", unique(x$plausible_range), "\n")
  if ("group" %in% colnames(x)) {
    groups <- unique(x$group)
    for (g in groups) {
      sub <- x[x$group == g, ]
      pct <- round(100 * mean(sub$informative), 1)
      cat(sprintf("  [%s]  %.1f%% informative  mean ratio = %.3f\n",
                  g, pct, mean(sub$ratio)))
    }
  } else {
    cat("  Informative      :", sum(x$informative), "/", nrow(x), "\n")
    cat("  Mean CI/range    :", round(mean(x$ratio), 3), "\n")
    cat("  Max CI/range     :", round(max(x$ratio), 3), "\n")
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
