# R/plot.R — et_plot_*(): ggplot2 visualizations (all return ggplot objects)

# ******************************************************************************
# et_plot_decomposition(): stacked bar of variance components
# ______________________________________________________________________________

#' Plot uncertainty decomposition
#'
#' Produces a stacked bar chart showing the relative contributions of
#' parameter, environmental, and residual variance for each observation,
#' plus a fourth temporal-autocorrelation component when present in
#' \code{decomp} (see \code{\link{decompose_uncertainty}}).
#'
#' @param decomp A \code{data.frame} from \code{\link{decompose_uncertainty}}
#'   or directly the \code{$decomposition} slot of an \code{et_prediction}.
#' @param proportional Logical.  If \code{TRUE} (default), each component is
#'   shown as a share of \code{total_var} (the measured predictive variance).
#'   Shares therefore sum to 1 exactly for autocorrelation models and to within
#'   Monte Carlo error for iid ones, where the shortfall is the
#'   \code{closure_gap} column.  If \code{FALSE}, raw variances are shown.
#'   Because \code{temporal_var} is signed, a negative component is drawn below
#'   the zero line rather than silently clamped.
#' @param group_col Character.  Optional name of a grouping column in
#'   \code{decomp} (present for grouped predictions).  One bar is drawn per
#'   group: the components are summed over the group's observations and
#'   divided by the summed \code{total_var}, so the bar still sums to 1.  That
#'   is a variance-weighted share, equal to
#'   \code{mean(component) / mean(total_var)}, and it is the quantity to quote
#'   alongside the figure.  It is not the mean of the per-observation shares,
#'   which would weight a low-variance lead equally with a high-variance one;
#'   for an autocorrelation model, where predictive variance grows with lead
#'   time, the two can differ by tens of percentage points.
#' @return A \code{ggplot2} object.
#' @export
et_plot_decomposition <- function(decomp, proportional = TRUE,
                                   group_col = NULL) {
  components <- c("param_var",
                  if ("group_var" %in% colnames(decomp)) "group_var",
                  "env_var", "residual_var",
                  if ("temporal_var" %in% colnames(decomp)) "temporal_var")

  grouped <- !is.null(group_col) && group_col %in% colnames(decomp)

  # One bar per group puts many observations on one x position. Stacking their
  # per-observation shares would make the bar sum to the NUMBER OF OBSERVATIONS
  # rather than to 1 -- the y axis then reads 16 on a 16-step window while
  # still being labelled a share. Aggregating first is also the only way the
  # bar can match the share reported in prose, which divides the mean component
  # by the mean total rather than averaging the per-observation ratios.
  if (grouped) {
    keep <- c(components, "total_var")
    keep <- keep[keep %in% colnames(decomp)]
    agg  <- stats::aggregate(decomp[keep],
                             by  = list(decomp[[group_col]]),
                             FUN = function(x) sum(x, na.rm = TRUE))
    names(agg)[1]    <- group_col
    agg[["obs_id"]]  <- seq_len(nrow(agg))
    decomp           <- agg
  }

  # Normalise against the MEASURED total rather than against the component sum.
  # position = "fill" would divide by the sum of the plotted components, which
  # (a) hides the iid closure gap by rescaling it away and (b) is meaningless
  # once a component can be negative. Dividing by total_var keeps the y axis
  # interpretable as "share of predictive variance" in every case.
  if (proportional && "total_var" %in% colnames(decomp)) {
    denom <- decomp[["total_var"]]
    denom[!is.finite(denom) | denom <= 0] <- NA_real_
    for (cc in components) decomp[[cc]] <- decomp[[cc]] / denom
  }

  long <- .pivot_decomp(decomp, components)

  has_negative <- any(long[["value"]] < 0, na.rm = TRUE)
  if (has_negative) {
    .et_info("A variance component is negative (temporal_var is a signed ",
             "covariance term); it is drawn below the zero line. A negative ",
             "temporal_var of the same order as v_env_mcse is Monte Carlo ",
             "noise rather than a real reduction in predictive variance.")
  }

  x_col <- if (!is.null(group_col) && group_col %in% colnames(decomp)) {
    group_col
  } else {
    "obs_id"
  }

  fill_values <- c(
    "Parameter"     = "#2166AC",
    "Environmental" = "#4DAC26",
    "Residual"      = "#D1E5F0",
    "Temporal"      = "#B2182B",
    "Group"         = "#762A83"
  )

  p <- ggplot2::ggplot(long, ggplot2::aes(
    x    = .data[[x_col]],
    y    = .data[["value"]],
    fill = .data[["component"]]
  )) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = fill_values, name = "Component") +
    ggplot2::labs(
      title = "Uncertainty Decomposition",
      x     = if (x_col == "obs_id") "Observation" else x_col,
      y     = if (proportional) "Share of predictive variance" else "Variance"
    ) +
    et_theme()

  # Angled ticks earn their place for a long run of observation ids; for a
  # handful of group names they are just harder to read.
  if (x_col == "obs_id") {
    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }

  # Make the zero line visible whenever a signed component can cross it.
  if (has_negative) {
    p <- p + ggplot2::geom_hline(yintercept = 0, linewidth = 0.3,
                                 colour = "grey20")
  }

  p
}

# ______________________________________________________________________________
# Where to put the threshold tag.
#
# A fixed anchor at the left edge is wrong as often as it is right: on a panel
# carrying several series the left edge is usually where they start, so the tag
# lands on top of them and the reviewer cannot read it. This scores a grid of x
# positions by how much clear vertical room there is beside the threshold line,
# above it and below it, counting only points that are actually VISIBLE (a
# capped panel hides the rest), and returns the roomiest. Ties keep the
# leftmost position, which is where a reader looks first.
.threshold_anchor <- function(x, y, hline, ylim = NULL, n_grid = 25L) {
  fallback <- list(where = "panel",
                   x = if (length(x)) min(x, na.rm = TRUE) else 0,
                   hjust = 0, vjust = -0.2)
  ok <- is.finite(x) & is.finite(y)
  if (!is.null(ylim)) ok <- ok & y >= min(ylim) & y <= max(ylim)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L) return(fallback)
  rng <- range(x)
  if (!is.finite(diff(rng)) || diff(rng) <= 0) return(fallback)

  # Room to the panel edge when nothing lies on that side at all.
  top <- if (!is.null(ylim)) max(ylim) else max(y, hline)
  bot <- if (!is.null(ylim)) min(ylim) else min(y, hline)

  grid <- seq(rng[1L], rng[2L], length.out = n_grid)
  win  <- diff(rng) * 0.10
  best <- list(score = -Inf)
  for (gx in grid) {
    near <- abs(x - gx) <= win
    if (!any(near)) next
    yy    <- y[near]
    above <- yy[yy > hline]
    below <- yy[yy < hline]
    room_above <- if (length(above)) min(above) - hline else top - hline
    room_below <- if (length(below)) hline - max(below) else hline - bot
    for (cand in list(list(score = room_above, vjust = -0.2),
                      list(score = room_below, vjust =  1.2))) {
      if (cand$score > best$score) {
        best <- list(score = cand$score, x = gx, vjust = cand$vjust)
      }
    }
  }
  if (!is.finite(best$score)) return(fallback)

  # Clearance is only useful measured against the size of the thing being
  # placed. The tag is about a tenth of the panel high, so a gap narrower than
  # that is not a gap: the series would run through the box. When no position
  # on the line clears it -- which is the normal case for several series that
  # all hover near the threshold, since that is what makes the threshold worth
  # drawing -- the tag goes outside the panel at the end of the line, where it
  # cannot collide with anything.
  y_span <- (if (!is.null(ylim)) diff(range(ylim)) else diff(range(c(y, hline))))
  if (!is.finite(y_span) || y_span <= 0) return(fallback)
  if (best$score < 0.10 * y_span) {
    return(list(where = "margin", x = rng[2L], hjust = 0, vjust = 0.5))
  }

  # Keep the box inside the panel: anchor it by whichever edge it is nearest.
  pos <- (best$x - rng[1L]) / diff(rng)
  best$hjust <- if (pos < 0.25) 0 else if (pos > 0.75) 1 else 0.5
  best$where <- "panel"
  best
}

# ******************************************************************************
# et_plot_shelf_life(): CI width vs time
# ______________________________________________________________________________

#' Plot forecast shelf life
#'
#' Shows how credible interval width grows over the forecast horizon and
#' marks the threshold beyond which the forecast is uninformative.
#'
#' @param sl An \code{et_shelf_life} object from \code{\link{shelf_life}}.
#' @param show_ratio Logical.  If \code{TRUE} (default), plots the
#'   ratio (CI width / plausible range) rather than raw CI width.
#' @param max_ratio Numeric.  Optional upper limit for the y axis, applied as a
#'   zoom (\code{coord_cartesian}) rather than a filter, so nothing is dropped
#'   and a trace that exceeds it runs off the top of the panel.  Use it when
#'   one series sits an order of magnitude above the rest and compresses them
#'   against the axis; say in the caption that the panel is capped.
#'   \code{NULL} (default) fits all of the data.
#' @return A \code{ggplot2} object.
#' @export
et_plot_shelf_life <- function(sl, show_ratio = TRUE, max_ratio = NULL) {
  y_col  <- if (show_ratio) "ratio"    else "ci_width"
  hline  <- if (show_ratio) 1.0        else unique(sl$plausible_range)[1]
  y_lab  <- if (show_ratio) "CI width / plausible range" else "CI width"
  hline_lab <- if (show_ratio) "Uninformative threshold"
               else "Plausible range"

  has_group <- "group" %in% colnames(sl)

  tag <- .threshold_anchor(sl[["time"]], sl[[y_col]], hline,
                           ylim = if (!is.null(max_ratio) && is.finite(max_ratio))
                                    c(0, max_ratio) else NULL)

  p <- ggplot2::ggplot(sl, ggplot2::aes(
    x     = .data[["time"]],
    y     = .data[[y_col]],
    color = if (has_group) .data[["group"]] else NULL,
    group = if (has_group) .data[["group"]] else NULL
  )) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_hline(yintercept = hline, linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      title = "Forecast Shelf Life",
      x     = "Time",
      y     = y_lab,
      color = if (has_group) "Group" else NULL
    ) +
    et_theme()

  if (tag$where == "panel") {
    p <- p + ggplot2::annotate("label", x = tag$x, y = hline,
                               label = hline_lab,
                               hjust = tag$hjust, vjust = tag$vjust,
                               colour = "grey20", fill = "white", size = 3.5)
  } else {
    # Outside the panel, as a break on a secondary axis at the line's own
    # value. Annotating into the margin would need clip = "off", which stops
    # clipping the data too and lets a steep series run across the title.
    p <- p +
      ggplot2::scale_y_continuous(
        sec.axis = ggplot2::sec_axis(~ ., breaks = hline,
                                     labels = sub(" ", "\n", hline_lab))) +
      ggplot2::theme(
        axis.text.y.right  = ggplot2::element_text(colour = "grey20", size = 8,
                                                   lineheight = 0.95, hjust = 0),
        axis.ticks.y.right = ggplot2::element_blank(),
        axis.title.y.right = ggplot2::element_blank())
  }

  if (has_group) p <- p + ggplot2::scale_color_viridis_d()

  # A zoom, not a filter: the traces still run to the panel edge.
  if (!is.null(max_ratio) && is.finite(max_ratio)) {
    p <- p + ggplot2::coord_cartesian(ylim = c(0, max_ratio))
  }

  p
}

# ******************************************************************************
# et_plot_calibration(): observed vs nominal coverage
# ______________________________________________________________________________

#' Plot calibration: observed vs nominal coverage
#'
#' A well-calibrated model produces points along the 1:1 diagonal.
#' Points above the diagonal indicate over-coverage (conservative);
#' below indicates under-coverage (anti-conservative).
#'
#' @param cal A \code{data.frame} from \code{\link{et_calibrate}}.
#' @param group_col Optional character.  Name of a column in \code{cal} that
#'   identifies sub-groups (e.g.\ \code{"species"}, \code{"cluster_id"}).
#'   When \code{NULL} (default), \code{et_plot_calibration} uses the
#'   \code{group} column if present; otherwise it auto-detects any single
#'   non-canonical column with more than one unique value and treats it as
#'   the grouping.  Set to \code{NA} to force a single un-grouped series.
#' @return A \code{ggplot2} object.
#' @export
et_plot_calibration <- function(cal, group_col = NULL) {
  canonical <- c("ci_level", "nominal", "observed_coverage",
                 "n_obs", "calibration_error", "sharpness")

  if (is.null(group_col)) {
    if ("group" %in% colnames(cal)) {
      group_col <- "group"
    } else {
      extras <- setdiff(colnames(cal), canonical)
      extras <- extras[vapply(extras,
                              function(x) length(unique(cal[[x]])) > 1,
                              logical(1))]
      if (length(extras) == 1L) group_col <- extras
    }
  } else if (is.na(group_col)) {
    group_col <- NULL
  }

  has_group <- !is.null(group_col) && group_col %in% colnames(cal)

  p <- ggplot2::ggplot(cal, ggplot2::aes(
    x = .data[["nominal"]],
    y = .data[["observed_coverage"]]
  )) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", color = "grey50") +
    ggplot2::labs(
      title = "Calibration: Observed vs Nominal Coverage",
      x     = "Nominal CI level",
      y     = "Observed coverage"
    ) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    et_theme()

  if (has_group) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(
        color = .data[[group_col]], group = .data[[group_col]]
      ), alpha = 0.6) +
      ggplot2::geom_point(ggplot2::aes(color = .data[[group_col]]), size = 3) +
      ggplot2::scale_color_viridis_d(name = group_col)
  } else {
    p <- p +
      ggplot2::geom_line(alpha = 0.7) +
      ggplot2::geom_point(size = 3, color = "steelblue")
  }

  p
}

# ******************************************************************************
# et_plot_forecast(): fan chart of posterior predictive intervals
# ______________________________________________________________________________

#' Plot posterior predictive fan chart
#'
#' Overlays nested credible interval ribbons on the median forecast, with
#' optional observed values for calibration assessment.
#'
#' @param predictions An \code{et_prediction} object.
#' @param observed Optional \code{data.frame} with true response values.
#'   If provided, points are overlaid.
#' @param response_col Character.  Name of the response column in
#'   \code{observed}.
#' @param time_col Character.  Column in \code{predictions$newdata} used
#'   as the x-axis.  Defaults to observation index.
#' @return A \code{ggplot2} object.
#' @export
et_plot_forecast <- function(predictions, observed = NULL,
                              response_col = NULL, time_col = NULL) {
  ci_df <- predictions$credible_intervals

  time_vals <- if (!is.null(time_col) &&
                   time_col %in% colnames(predictions$newdata)) {
    rep(predictions$newdata[[time_col]],
        times = length(unique(ci_df$ci_level)))
  } else {
    ci_df$row_id
  }
  ci_df$time <- time_vals

  p <- ggplot2::ggplot(ci_df, ggplot2::aes(x = .data[["time"]])) +
    ggplot2::geom_ribbon(ggplot2::aes(
      ymin  = .data[["lower"]],
      ymax  = .data[["upper"]],
      alpha = .data[["ci_level"]],
      group = .data[["ci_level"]]
    ), fill = "steelblue") +
    ggplot2::geom_line(
      data = ci_df[ci_df$ci_level == max(ci_df$ci_level), ],
      ggplot2::aes(y = .data[["median"]]),
      color = "darkblue", linewidth = 0.9
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.12, 0.45), guide = "none") +
    ggplot2::labs(
      title = "Forecast Fan Chart",
      x     = if (!is.null(time_col)) time_col else "Observation",
      y     = "Response"
    ) +
    et_theme()

  if (!is.null(observed) && !is.null(response_col)) {
    obs_time <- if (!is.null(time_col) &&
                    time_col %in% colnames(observed)) {
      observed[[time_col]]
    } else {
      seq_len(nrow(observed))
    }
    obs_df <- data.frame(
      time     = obs_time,
      response = observed[[response_col]],
      stringsAsFactors = FALSE
    )
    p <- p + ggplot2::geom_point(
      data = obs_df,
      ggplot2::aes(x = .data[["time"]], y = .data[["response"]]),
      color = "firebrick", size = 2.5, inherit.aes = FALSE
    )
  }

  p
}

# ******************************************************************************
# et_plot_prior_posterior(): density overlay
# ______________________________________________________________________________

#' Plot prior vs posterior distributions for model coefficients
#'
#' Overlays prior and posterior density for each predictor coefficient,
#' visualising how much the data update the priors.
#'
#' @param model An \code{et_model} object.
#' @param max_preds Integer.  Maximum number of predictors to show
#'   (default 8).  Predictors are shown in the order they appear in the
#'   prior specification.
#' @param n_prior_draws Integer.  Number of random draws for the prior
#'   density (default 4000).
#' @return A \code{ggplot2} object.
#' @export
et_plot_prior_posterior <- function(model, max_preds = 8L,
                                     n_prior_draws = 4000L) {
  if (!inherits(model, "et_model")) {
    stop("model must be an et_model object.")
  }
  if (is.null(model$prior_spec)) {
    stop("No et_prior_spec attached to model. ",
         "Prior/posterior plot requires a model fitted via et_fit() with ",
         "priors from extract_priors().")
  }

  fit     <- model$fit
  spec    <- model$prior_spec
  preds   <- head(spec$pred_names, max_preds)

  post_draws <- as.data.frame(fit)
  beta_cols  <- paste0("b_", preds)
  avail      <- intersect(beta_cols, colnames(post_draws))
  if (length(avail) == 0) stop("No beta columns found in posterior draws.")

  post_long <- do.call(rbind, lapply(avail, function(bc) {
    pname <- sub("^b_", "", bc)
    data.frame(
      parameter    = pname,
      value        = post_draws[[bc]],
      distribution = "Posterior",
      stringsAsFactors = FALSE
    )
  }))

  prior_long <- do.call(rbind, lapply(avail, function(bc) {
    pname <- sub("^b_", "", bc)
    idx   <- match(pname, spec$pred_names)
    mu    <- if (!is.null(spec$coefs)) spec$coefs[idx] else 0
    sd    <- max(spec$multiplier * abs(mu), spec$min_sd)
    data.frame(
      parameter    = pname,
      value        = stats::rnorm(n_prior_draws, mu, sd),
      distribution = "Prior",
      stringsAsFactors = FALSE
    )
  }))

  combined <- rbind(post_long, prior_long)
  combined$distribution <- factor(combined$distribution, levels = c("Prior", "Posterior"))

  ggplot2::ggplot(combined,
                  ggplot2::aes(x = .data[["value"]], fill = .data[["distribution"]])) +
    ggplot2::geom_density(alpha = 0.45) +
    ggplot2::facet_wrap(~ .data[["parameter"]], scales = "free") +
    ggplot2::scale_fill_manual(
      values = c("Prior" = "grey60", "Posterior" = "steelblue"),
      name   = NULL
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "grey30") +
    ggplot2::labs(
      title = "Prior vs Posterior: Regression Coefficients",
      x     = "Value",
      y     = "Density"
    ) +
    et_theme()
}

# ******************************************************************************
# et_plot_coefficients(): forest plot, Bayesian vs regularized
# ______________________________________________________________________________

#' Forest plot of regression coefficients
#'
#' Compares Bayesian posterior estimates (95\% CI) with the regularized
#' coefficient values used as prior means (shown as crosses).
#'
#' @param model An \code{et_model} or \code{et_model_list}.
#' @return A \code{ggplot2} object.
#' @export
et_plot_coefficients <- function(model) {
  UseMethod("et_plot_coefficients")
}

#' @export
et_plot_coefficients.et_model <- function(model) {
  coef_df <- .extract_coef_df(model, group_label = NULL)
  .forest_plot(coef_df, facet_by_group = FALSE)
}

#' @export
et_plot_coefficients.et_model_list <- function(model) {
  parts <- lapply(names(model$models), function(g) {
    m <- model$models[[g]]
    if (is.null(m)) return(NULL)
    .extract_coef_df(m, group_label = g)
  })
  coef_df <- do.call(rbind, Filter(Negate(is.null), parts))
  .forest_plot(coef_df, facet_by_group = TRUE)
}

.extract_coef_df <- function(model, group_label) {
  fit  <- model$fit
  spec <- model$prior_spec
  post <- brms::fixef(fit)
  fe   <- rownames(post)
  fe   <- fe[fe != "Intercept"]

  rows <- lapply(fe, function(pname) {
    enet_val <- if (!is.null(spec) && !is.null(spec$coefs) &&
                    pname %in% names(spec$coefs)) {
      spec$coefs[pname]
    } else {
      NA_real_
    }
    df <- data.frame(
      predictor  = pname,
      post_mean  = post[pname, "Estimate"],
      post_lower = post[pname, "Q2.5"],
      post_upper = post[pname, "Q97.5"],
      enet_coef  = enet_val,
      stringsAsFactors = FALSE
    )
    if (!is.null(group_label)) df$group <- group_label
    df
  })
  do.call(rbind, rows)
}

.forest_plot <- function(coef_df, facet_by_group) {
  p <- ggplot2::ggplot(coef_df,
    ggplot2::aes(y = .data[["predictor"]])) +
    ggplot2::geom_pointrange(
      ggplot2::aes(
        x    = .data[["post_mean"]],
        xmin = .data[["post_lower"]],
        xmax = .data[["post_upper"]]
      ),
      color = "steelblue", linewidth = 0.6
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    ggplot2::labs(
      title = "Coefficient Estimates: Posterior (blue) vs Regularized (red \u00d7)",
      x     = "Coefficient",
      y     = NULL
    ) +
    et_theme()

  # Add enet estimates if present
  has_enet <- any(!is.na(coef_df$enet_coef))
  if (has_enet) {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(x = .data[["enet_coef"]]),
      color = "firebrick", shape = 4, size = 3, na.rm = TRUE
    )
  }

  if (facet_by_group && "group" %in% colnames(coef_df)) {
    p <- p + ggplot2::facet_wrap(~ .data[["group"]], scales = "free")
  }

  p
}

# ******************************************************************************
# et_plot_sensitivity(): noise fraction vs horizon / env_share
# ______________________________________________________________________________

#' Plot a sensitivity profile
#'
#' Visualises the output of \code{\link{et_sensitivity_profile}}: for each
#' noise grid point, shows how the \strong{environmental share} of total
#' variance and the \strong{forecast horizon} respond.  Observed horizons
#' are drawn as solid points; projected horizons are hollow points;
#' lower-bound rows are drawn as upward arrows at the last informative time.
#'
#' The x-axis is the noise fraction (when the grid was built from
#' \code{fraction_grid}) or the grid step label otherwise.  When both a
#' numeric fraction and a descriptive label exist, the fraction is preferred
#' for continuous x-positioning.
#'
#' @param sens An \code{et_sensitivity} object from
#'   \code{\link{et_sensitivity_profile}}.
#' @param show \code{"horizon"} (default) shows the shelf-life horizon;
#'   \code{"env_share"} shows \code{env_var / total_var}; \code{"ratio"}
#'   shows the mean CI width / plausible range.
#' @return A \code{ggplot2} object.
#' @seealso \code{\link{et_sensitivity_profile}}
#' @export
et_plot_sensitivity <- function(sens, show = c("horizon", "env_share", "ratio")) {
  if (!inherits(sens, "et_sensitivity")) {
    stop("et_plot_sensitivity() expects an et_sensitivity object.")
  }
  show <- match.arg(show)

  df <- as.data.frame(sens)
  x_is_numeric <- !all(is.na(df$fraction))
  df$x <- if (x_is_numeric) df$fraction else seq_len(nrow(df))
  x_lab <- if (x_is_numeric) "Noise fraction of predictor SD" else "Grid step"

  if (show == "horizon") {
    df$y <- df$horizon
    df$point_shape <- factor(
      ifelse(df$horizon_type == "observed", "observed",
      ifelse(df$horizon_type == "projected", "projected", "lower bound")),
      levels = c("observed", "projected", "lower bound")
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]])) +
      ggplot2::geom_line(na.rm = TRUE, colour = "#2166AC") +
      ggplot2::geom_point(ggplot2::aes(shape = .data[["point_shape"]]),
                          size = 3, na.rm = TRUE) +
      ggplot2::scale_shape_manual(
        values = c("observed" = 16, "projected" = 21, "lower bound" = 2),
        drop = FALSE, name = "Horizon type"
      ) +
      ggplot2::labs(x = x_lab, y = "Forecast horizon",
                    title = "Shelf life vs. assumed measurement noise")
  } else if (show == "env_share") {
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["x"]],
                                          y = .data[["env_share"]])) +
      ggplot2::geom_line(colour = "#4DAC26") +
      ggplot2::geom_point(colour = "#4DAC26", size = 3) +
      ggplot2::scale_y_continuous(labels = scales_percent_safe,
                                  limits = c(0, NA)) +
      ggplot2::labs(x = x_lab,
                    y = "Environmental share of total variance",
                    title = "How much of the forecast variance is environmental?")
  } else {
    thr <- attr(sens, "threshold")
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["x"]],
                                          y = .data[["ratio_mean"]])) +
      ggplot2::geom_line(colour = "#B2182B") +
      ggplot2::geom_point(colour = "#B2182B", size = 3) +
      ggplot2::geom_hline(yintercept = if (is.null(thr)) 1 else thr,
                          linetype = "dashed") +
      ggplot2::labs(x = x_lab,
                    y = "Mean CI width / plausible range",
                    title = "Forecast informativeness vs. noise")
  }

  p + et_theme()
}

#' Plot a PIT / rank histogram for calibration
#'
#' Histogram of the probability-integral-transform values from
#' \code{\link{et_pit}}.  Bars level with the dashed uniform reference line
#' indicate calibration; a \strong{U-shape} means overconfident (too-narrow)
#' intervals, a central \strong{hump} means underconfident (too-wide), and a
#' \strong{tilt} means bias.  The shaded band is an approximate 95\% consistency
#' region for the bar heights under uniformity.  When \code{pit} carries a
#' \code{group} column the reference line and band are computed separately for
#' each panel, so panels holding different numbers of observations remain
#' comparable.
#'
#' @param pit A \code{data.frame} from \code{\link{et_pit}} (a \code{pit} column,
#'   optionally \code{group}), or a bare numeric vector of PIT values.
#' @param bins Integer.  Number of histogram bins (default 10).
#' @return A \code{ggplot} object.
#' @seealso \code{\link{et_pit}}, \code{\link{et_plot_calibration}}
#' @export
et_plot_pit <- function(pit, bins = 10L) {
  df <- if (is.numeric(pit)) data.frame(pit = pit) else pit
  if (!"pit" %in% colnames(df)) {
    stop("`pit` must be a numeric vector or a data.frame with a 'pit' column ",
         "(as returned by et_pit()).")
  }
  bins      <- max(2L, as.integer(bins))
  has_group <- "group" %in% colnames(df)
  breaks    <- seq(0, 1, length.out = bins + 1L)

  # Expected bar height and an approximate 95% band under uniformity, one per
  # panel. Groups routinely differ in size -- a 35-year forecast window beside
  # a 16-year one -- and a single pooled reference lands below the large
  # panel's bars and above the small panels', so a panel whose every bar
  # clears the line reads as grossly miscalibrated when a histogram cannot be:
  # its bars must average to the panel's own n / bins.
  ref <- function(n) {
    e <- n / bins
    b <- 1.96 * sqrt(e * (1 - 1 / bins))
    data.frame(expected = e, ymin = pmax(0, e - b), ymax = e + b, pit = 0)
  }
  if (has_group) {
    n_per <- tapply(!is.na(df$pit), df$group, sum)
    n_per <- n_per[!is.na(n_per) & n_per > 0]
    bands <- ref(as.numeric(n_per))
    bands$group <- if (is.factor(df$group)) {
      factor(names(n_per), levels = levels(df$group))
    } else {
      names(n_per)
    }
  } else {
    bands <- ref(sum(!is.na(df$pit)))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["pit"]])) +
    ggplot2::geom_rect(data = bands, inherit.aes = FALSE,
                       ggplot2::aes(ymin = .data[["ymin"]],
                                    ymax = .data[["ymax"]]),
                       xmin = 0, xmax = 1, fill = "grey80", alpha = 0.5) +
    ggplot2::geom_histogram(breaks = breaks, boundary = 0,
                            fill = "steelblue", color = "white") +
    ggplot2::geom_hline(data = bands, inherit.aes = FALSE,
                        ggplot2::aes(yintercept = .data[["expected"]]),
                        linetype = "dashed", color = "grey40") +
    ggplot2::labs(x = "PIT", y = "Count",
                  title = "PIT histogram (calibration)") +
    et_theme()

  if (has_group) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[["group"]]))
  p
}

# Permissive wrapper for scales::percent — fall back to identity if the
# scales package isn't installed so the plot still renders.
scales_percent_safe <- function(x) {
  if (requireNamespace("scales", quietly = TRUE)) scales::percent(x)
  else formatC(100 * x, format = "f", digits = 0, flag = "#")
}

# ******************************************************************************
# Internal helpers
# ______________________________________________________________________________

# Pivot variance decomposition to long format for ggplot stacking
.pivot_decomp <- function(decomp, cols) {
  rows <- lapply(cols, function(col) {
    label <- switch(col,
      param_var    = "Parameter",
      env_var      = "Environmental",
      residual_var = "Residual",
      temporal_var = "Temporal",
      group_var    = "Group",
      col
    )
    df <- decomp[, c(setdiff(colnames(decomp), cols), col), drop = FALSE]
    df$component <- label
    df$value     <- decomp[[col]]
    df[, c(setdiff(colnames(decomp), cols), "component", "value")]
  })
  do.call(rbind, rows)
}
