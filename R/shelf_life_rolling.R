# R/shelf_life_rolling.R -- pooling shelf lives across forecast origins

#' Pool forecast shelf lives across multiple forecast origins
#'
#' A shelf life read off a \emph{single} forecast origin is a property of that
#' training window as much as of the system: it inherits whatever was
#' idiosyncratic about the final training years.  Established practice in
#' ecological forecast evaluation is to issue forecasts repeatedly from rolling
#' origins and summarise skill indexed by \strong{both} origin and lead time
#' (Dietze 2017; Thomas et al. 2020; Petchey et al. 2015).  This function does
#' the pooling half of that: given the horizons obtained at several origins, it
#' returns a horizon with an interval rather than a point.
#'
#' @section Why censoring matters:
#' Later origins have fewer lead times available, so a long horizon may not be
#' observed before the data run out --- the forecast is still informative at the
#' last available lead.  That is \strong{right-censoring}, not a missing value.
#' Dropping those origins would bias the pooled horizon downward (only the short
#' horizons survive), and treating the last lead as the horizon would bias it
#' downward too.  Both are avoided by treating the pooled horizon as a
#' time-to-event problem and estimating it with Kaplan--Meier, which uses the
#' censored origins as the partial information they are.
#'
#' When \pkg{survival} is available the median horizon and its confidence
#' interval come from \code{survival::survfit()}.  Otherwise the function falls
#' back to the median of the observed horizons and says so in the
#' \code{method} field --- a fallback that is biased low whenever anything is
#' censored, which is why the interval is reported as \code{NA} there.
#'
#' @param horizons Either a list of \code{et_shelf_life} objects (one per
#'   forecast origin, in any order), or a \code{data.frame} with columns
#'   \code{lead} (numeric horizon in lead time) and \code{censored} (logical,
#'   \code{TRUE} when the forecast was still informative at the last available
#'   lead).  A \code{group} column, when present, is summarised separately.
#' @param conf_level Confidence level for the Kaplan--Meier interval on the
#'   median horizon (default \code{0.95}).
#'
#' @return A \code{data.frame}, one row per group (or a single row when
#'   ungrouped), with columns:
#'   \describe{
#'     \item{group}{Group label, when supplied.}
#'     \item{n_origins}{Number of forecast origins contributing.}
#'     \item{n_observed}{Origins at which the horizon was reached in-window.}
#'     \item{n_censored}{Origins right-censored (still informative at the end).}
#'     \item{median_lead}{Pooled horizon in lead time. \code{NA} when every
#'       origin is censored --- which is itself a result: the forecast never
#'       became uninformative from any origin.}
#'     \item{lower, upper}{Confidence bounds on the median.}
#'     \item{min_lead, max_lead}{Range of the observed horizons.}
#'     \item{method}{\code{"kaplan-meier"} or \code{"median (uncensored only)"}.}
#'   }
#'
#' @examples
#' # Horizons from five origins; the last was censored.
#' h <- data.frame(lead     = c(2, 3, 2, 3, 5),
#'                 censored = c(FALSE, FALSE, FALSE, FALSE, TRUE))
#' et_shelf_life_pool(h)
#'
#' @seealso \code{\link{shelf_life}} for the single-origin horizon.
#' @references
#' Dietze, M. C. (2017) Prediction in ecology: a first-principles framework.
#' \emph{Ecological Applications} 27:2048--2060.
#'
#' Petchey, O. L. \emph{et al.} (2015) The ecological forecast horizon, and
#' examples of its uses and determinants. \emph{Ecology Letters} 18:597--611.
#' @export
et_shelf_life_pool <- function(horizons, conf_level = 0.95) {
  df <- .as_horizon_frame(horizons)

  if (!nrow(df)) {
    stop("No forecast origins supplied.")
  }
  if (!all(c("lead", "censored") %in% colnames(df))) {
    stop("`horizons` must provide `lead` and `censored` columns.")
  }

  groups <- if ("group" %in% colnames(df)) as.character(df$group) else
    rep("__all__", nrow(df))

  out <- do.call(rbind, lapply(split(df, groups), function(d) {
    .pool_one_group(d, conf_level = conf_level)
  }))
  rownames(out) <- NULL

  if (!"group" %in% colnames(df)) out$group <- NULL
  out
}

# Coerce the accepted input forms to a data.frame of (lead, censored[, group]).
.as_horizon_frame <- function(horizons) {
  if (is.data.frame(horizons)) return(horizons)

  if (is.list(horizons)) {
    rows <- lapply(seq_along(horizons), function(i) {
      sl <- horizons[[i]]
      hor <- attr(sl, "horizon")
      if (is.null(hor)) {
        stop("Element ", i, " of `horizons` has no horizon attribute; ",
             "supply objects returned by shelf_life().")
      }
      observed <- identical(hor$type, "observed")
      # Lead time = position of the crossing among the forecast rows. Falls back
      # to the number of rows when censored (the last lead still informative).
      lead <- if (observed && !is.null(hor$value)) {
        match(hor$value, sort(unique(sl$time)))
      } else {
        length(unique(sl$time))
      }
      data.frame(origin   = names(horizons)[i] %||% as.character(i),
                 lead     = as.numeric(lead),
                 censored = !observed,
                 stringsAsFactors = FALSE)
    })
    return(do.call(rbind, rows))
  }

  stop("`horizons` must be a list of et_shelf_life objects or a data.frame.")
}

.pool_one_group <- function(d, conf_level) {
  obs   <- d[!d$censored, , drop = FALSE]
  n_obs <- nrow(obs)

  base <- data.frame(
    group      = d$group[1] %||% NA_character_,
    n_origins  = nrow(d),
    n_observed = n_obs,
    n_censored = sum(d$censored),
    median_lead = NA_real_,
    lower       = NA_real_,
    upper       = NA_real_,
    min_lead    = if (n_obs) min(obs$lead) else NA_real_,
    max_lead    = if (n_obs) max(obs$lead) else NA_real_,
    method      = NA_character_,
    stringsAsFactors = FALSE
  )
  if (!"group" %in% colnames(d)) base$group <- "__all__"

  # Every origin censored: the forecast never became uninformative. That is a
  # result, not a failure to estimate, and Kaplan-Meier has no median for it.
  if (n_obs == 0L) {
    base$method <- "all origins censored (no horizon reached)"
    return(base)
  }

  if (requireNamespace("survival", quietly = TRUE)) {
    km <- tryCatch({
      su <- survival::Surv(time = d$lead, event = !d$censored)
      sf <- survival::survfit(su ~ 1, conf.int = conf_level)
      tb <- summary(sf)$table
      # survfit names the bounds "<conf.int>LCL"/"<conf.int>UCL" on the
      # PROBABILITY scale, e.g. "0.95LCL" -- not "95LCL".
      lcl <- grep("LCL$", names(tb), value = TRUE)[1]
      ucl <- grep("UCL$", names(tb), value = TRUE)[1]
      list(med = unname(tb[["median"]]),
           lo  = if (!is.na(lcl)) unname(tb[[lcl]]) else NA_real_,
           hi  = if (!is.na(ucl)) unname(tb[[ucl]]) else NA_real_)
    }, error = function(e) NULL)

    if (!is.null(km)) {
      base$median_lead <- km$med
      base$lower       <- km$lo
      base$upper       <- km$hi
      base$method      <- "kaplan-meier"
      return(base)
    }
  }

  # Fallback: median of the observed horizons only. Biased low when anything is
  # censored, so the interval is deliberately left NA rather than computed from
  # a sample that excludes the long horizons.
  base$median_lead <- stats::median(obs$lead)
  base$method      <- if (any(d$censored)) {
    "median (uncensored only -- biased low; install 'survival' for Kaplan-Meier)"
  } else {
    "median (no censoring)"
  }
  base
}
