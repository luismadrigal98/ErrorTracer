# tests/testthat/test-plot.R
# Smoke tests for et_plot_* functions — verify they return ggplot objects
# and do not error with mock data, without requiring a real brms fit.

library(ggplot2)

# ── Helper: minimal mock decomposition data.frame ───────────────────────────

.mock_decomp <- function(n_obs = 5) {
  data.frame(
    obs_id       = seq_len(n_obs),
    total_var    = runif(n_obs, 0.3, 0.5),
    param_var    = runif(n_obs, 0.05, 0.15),
    env_var      = runif(n_obs, 0.02, 0.08),
    residual_var = rep(0.22, n_obs)
  )
}

.mock_shelf_life <- function(n_obs = 5) {
  sl <- data.frame(
    obs_id          = seq_len(n_obs),
    time            = 2020 + seq_len(n_obs) - 1L,
    ci_width        = runif(n_obs, 0.8, 1.4),
    plausible_range = 3.5,
    ratio           = runif(n_obs, 0.2, 0.4),
    informative     = TRUE,
    stringsAsFactors = FALSE
  )
  attr(sl, "threshold") <- 1.0
  attr(sl, "horizon")   <- list(
    value            = NA_real_,
    type             = "lower_bound",
    last_informative = 2024,
    description      = "All periods informative."
  )
  structure(sl, class = c("et_shelf_life", "data.frame"))
}

.mock_calibration <- function() {
  data.frame(
    ci_level           = c(0.50, 0.80, 0.90, 0.95),
    nominal            = c(0.50, 0.80, 0.90, 0.95),
    observed_coverage  = c(0.60, 0.80, 1.00, 1.00),
    n_obs              = 5L,
    calibration_error  = c(0.10, 0.00, 0.10, 0.05),
    stringsAsFactors   = FALSE
  )
}

# ── et_plot_decomposition ────────────────────────────────────────────────────

test_that("et_plot_decomposition returns a ggplot object (proportional)", {
  p <- et_plot_decomposition(.mock_decomp(), proportional = TRUE)
  expect_s3_class(p, "ggplot")
})

test_that("et_plot_decomposition returns a ggplot object (absolute)", {
  p <- et_plot_decomposition(.mock_decomp(), proportional = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("et_plot_decomposition handles grouped decomp with group_col", {
  d <- rbind(
    cbind(group = "A", .mock_decomp(3)),
    cbind(group = "B", .mock_decomp(3))
  )
  p <- et_plot_decomposition(d, group_col = "group")
  expect_s3_class(p, "ggplot")
})

# ── et_plot_shelf_life ───────────────────────────────────────────────────────

test_that("et_plot_shelf_life returns a ggplot object", {
  p <- et_plot_shelf_life(.mock_shelf_life())
  expect_s3_class(p, "ggplot")
})

test_that("et_plot_shelf_life with show_ratio = FALSE works", {
  p <- et_plot_shelf_life(.mock_shelf_life(), show_ratio = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("et_plot_shelf_life with grouped et_shelf_life adds colour", {
  sl_g <- cbind(group = rep(c("A","B"), each = 3),
                .mock_shelf_life(6))
  attr(sl_g, "threshold") <- 1.0
  class(sl_g) <- c("et_shelf_life", "data.frame")
  p <- et_plot_shelf_life(sl_g)
  expect_s3_class(p, "ggplot")
})

# ── et_plot_calibration ──────────────────────────────────────────────────────

test_that("et_plot_calibration returns a ggplot object (single group)", {
  p <- et_plot_calibration(.mock_calibration())
  expect_s3_class(p, "ggplot")
})

test_that("et_plot_calibration returns a ggplot object (grouped)", {
  cal <- rbind(
    cbind(group = "A", .mock_calibration()),
    cbind(group = "B", .mock_calibration())
  )
  p <- et_plot_calibration(cal)
  expect_s3_class(p, "ggplot")
})

# ── grouped decomposition: the bar is a share, not a sum of shares ───────────
# Regression test for a defect that survived the smoke test above: with one bar
# per group, geom_col stacked every observation's share on the same x position,
# so a 16-step window produced a bar of height 16 on an axis labelled "share of
# predictive variance". The axis was wrong; the colour proportions were not.

test_that("a grouped proportional bar sums to 1, not to the row count", {
  mk <- function(g, n = 16) {
    d <- .mock_decomp(n)
    d$total_var <- d$param_var + d$env_var + d$residual_var  # a closed budget
    cbind(group = g, d)
  }
  d <- rbind(mk("A"), mk("B"))
  p <- et_plot_decomposition(d, proportional = TRUE, group_col = "group")
  tot <- tapply(p$data$value, p$data$group, sum)
  expect_equal(as.numeric(tot), c(1, 1), tolerance = 1e-8)
  expect_equal(nrow(p$data), 2L * 3L)   # 2 groups x 3 components
})

test_that("a grouped bar uses the variance-weighted share quoted in prose", {
  # sum(component)/sum(total) == mean(component)/mean(total), and that is NOT
  # the mean of the per-observation ratios once total_var varies across the
  # window, which is exactly what an ar() model does.
  d <- data.frame(
    obs_id       = 1:3,
    group        = "A",
    total_var    = c(1, 10, 100),
    param_var    = c(0.5, 1, 5),
    env_var      = c(0.1, 0.5, 5),
    residual_var = c(0.4, 8.5, 90)
  )
  p <- et_plot_decomposition(d, proportional = TRUE, group_col = "group")
  got <- p$data$value[p$data$component == "Parameter"]
  expect_equal(got, mean(d$param_var) / mean(d$total_var), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(got, mean(d$param_var / d$total_var))))
})

test_that("ungrouped bars stay per-observation", {
  p <- et_plot_decomposition(.mock_decomp(5), proportional = TRUE)
  expect_equal(length(unique(p$data$obs_id)), 5L)
})

# ── et_plot_shelf_life(max_ratio): zoom, never drop ──────────────────────────

test_that("max_ratio zooms the panel without dropping data", {
  sl <- .mock_shelf_life(5)
  sl$ratio <- c(0.2, 0.5, 1.1, 3.0, 6.0)
  p <- et_plot_shelf_life(sl, max_ratio = 2)
  expect_equal(p$coordinates$limits$y, c(0, 2))
  expect_equal(nrow(p$data), 5L)                     # nothing filtered out
  expect_equal(max(p$data$ratio), 6.0)
})

# ── threshold tag placement ──────────────────────────────────────────────────
# The tag used to be pinned to the left edge at the threshold, which is exactly
# where a set of series that hover around that threshold starts, so it landed
# on the data. Placement is now chosen by clearance, and falls out of the panel
# altogether when nothing on the line clears the height of the tag.

test_that("the threshold tag stays on the line when there is room for it", {
  x <- 1:20
  y <- rep(0.2, 20)                     # everything far below the line
  a <- ErrorTracer:::.threshold_anchor(x, y, hline = 1, ylim = c(0, 2))
  expect_equal(a$where, "panel")
  expect_lt(a$vjust, 0)                 # sits above the line, where the room is
})

test_that("the threshold tag leaves the panel when the series hug the line", {
  x <- rep(1:20, 3)
  y <- c(rep(0.95, 20), rep(1.05, 20), rep(1.10, 20))
  a <- ErrorTracer:::.threshold_anchor(x, y, hline = 1, ylim = c(0, 2))
  expect_equal(a$where, "margin")
})

test_that("a hidden series does not count against the tag", {
  # A trace that leaves a capped panel cannot collide with anything.
  x <- rep(1:20, 2)
  y <- c(rep(0.2, 20), seq(3, 12, length.out = 20))
  a <- ErrorTracer:::.threshold_anchor(x, y, hline = 1, ylim = c(0, 2))
  expect_equal(a$where, "panel")
})

test_that("et_plot_shelf_life draws in both placements", {
  sl <- .mock_shelf_life(6)
  sl$ratio <- rep(0.2, 6)
  expect_s3_class(et_plot_shelf_life(sl), "ggplot")
  sl$ratio <- rep(1.02, 6)
  expect_s3_class(et_plot_shelf_life(sl), "ggplot")
})
