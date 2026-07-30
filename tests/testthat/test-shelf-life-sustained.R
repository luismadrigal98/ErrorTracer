# Regression tests for the sustained-crossing rule (ErrorTracer 1.3.1).
#
# Before 1.3.1 shelf_life() reported the FIRST period with ratio >= threshold as
# the horizon. When the ratio hovers near the threshold without trending, a
# single noisy period exceeds it by chance and was promoted to a "shelf life".
# This was observed on the Portal DM series, whose ratio sits at ~0.92 for the
# whole 16-year window with exactly one excursion to 1.038 -- reported as a
# 12-year horizon, and unstable to sampler settings alone.
#
# These tests operate directly on .compute_horizon(), so they need no Stan fit.

make_sl_df <- function(ratios, times = seq_along(ratios), threshold = 1.0) {
  data.frame(
    obs_id          = seq_along(ratios),
    time            = times,
    ci_width        = ratios,
    plausible_range = 1,
    ratio           = ratios,
    informative     = ratios < threshold
  )
}

test_that("an isolated excursion is not promoted to a horizon", {
  # The Portal DM hindcast trajectory: one value above 1.0 at index 12.
  dm <- c(0.935, 0.902, 0.894, 0.860, 0.921, 0.885, 0.965, 0.901,
          0.938, 0.975, 0.949, 1.038, 0.921, 0.854, 0.902, 0.942)
  h <- ErrorTracer:::.compute_horizon(make_sl_df(dm), threshold = 1.0,
                                      min_slope = 1e-4, min_run = 2L)
  expect_false(identical(h$type, "observed"))
  expect_equal(h$n_exceedances, 1L)
  expect_equal(h$first_exceedance, 12)
  # The excursion is still visible, just not treated as a crossing.
  expect_match(h$description, "isolated exceedance", fixed = TRUE)
})

test_that("min_run = 1 reproduces the pre-1.3.1 first-exceedance behaviour", {
  dm <- c(0.935, 0.902, 0.894, 0.860, 0.921, 0.885, 0.965, 0.901,
          0.938, 0.975, 0.949, 1.038, 0.921, 0.854, 0.902, 0.942)
  h <- ErrorTracer:::.compute_horizon(make_sl_df(dm), threshold = 1.0,
                                      min_slope = 1e-4, min_run = 1L)
  expect_identical(h$type, "observed")
  expect_equal(h$value, 12)
})

test_that("a genuine monotone degradation is still reported as observed", {
  # Portal PP: crosses at index 2 and never comes back.
  pp <- c(0.959, 1.133, 1.336, 1.566, 1.844, 2.113, 2.443, 2.667,
          3.043, 3.364, 3.724, 4.135, 4.516, 4.981, 5.564, 6.258)
  h <- ErrorTracer:::.compute_horizon(make_sl_df(pp), threshold = 1.0,
                                      min_slope = 1e-4, min_run = 2L)
  expect_identical(h$type, "observed")
  expect_equal(h$value, 2)
  expect_equal(h$last_informative, 1)
})

test_that("a crossing already in force at the first lead is caught", {
  # Portal DO: above threshold from lead 1 onward.
  do_r <- c(1.02, 1.09, 1.08, 1.14, 1.16, 1.19, 1.19, 1.12,
            1.17, 1.21, 1.12, 1.28, 1.19, 1.17, 1.17, 1.18)
  h <- ErrorTracer:::.compute_horizon(make_sl_df(do_r), threshold = 1.0,
                                      min_slope = 1e-4, min_run = 2L)
  expect_identical(h$type, "observed")
  expect_equal(h$value, 1)
  expect_true(is.na(h$last_informative))
})

test_that("a dip inside a degrading series does not reset the horizon", {
  # Portal DM reforecast: two early exceedances, one dip, then sustained.
  ref <- c(0.929, 0.950, 0.996, 1.014, 1.008, 0.998, 1.035, 1.033,
           1.055, 1.066, 1.078, 1.081, 1.099, 1.096, 1.146, 1.119)
  h2 <- ErrorTracer:::.compute_horizon(make_sl_df(ref), threshold = 1.0,
                                       min_slope = 1e-4, min_run = 2L)
  expect_identical(h2$type, "observed")
  expect_equal(h2$value, 4)          # the 1.014/1.008 pair persists for 2
  # A stricter run requirement skips the pair and lands on the sustained tail.
  h3 <- ErrorTracer:::.compute_horizon(make_sl_df(ref), threshold = 1.0,
                                       min_slope = 1e-4, min_run = 3L)
  expect_identical(h3$type, "observed")
  expect_equal(h3$value, 7)
})

test_that("min_run longer than the series yields no crossing", {
  h <- ErrorTracer:::.compute_horizon(make_sl_df(c(1.5, 1.6)), threshold = 1.0,
                                      min_slope = 1e-4, min_run = 5L)
  expect_false(identical(h$type, "observed"))
  expect_equal(h$n_exceedances, 2L)
})

test_that("diagnostics are attached to every horizon type", {
  flat <- rep(0.5, 10)                       # lower_bound, no exceedance
  h <- ErrorTracer:::.compute_horizon(make_sl_df(flat), threshold = 1.0,
                                      min_slope = 1e-4, min_run = 2L)
  expect_identical(h$type, "lower_bound")
  expect_equal(h$n_exceedances, 0L)
  expect_equal(h$frac_exceedance, 0)
  expect_true(is.na(h$first_exceedance))
  expect_equal(h$min_run, 2L)
})

test_that("a flat, noisy trend is not extrapolated into a projection", {
  # Portal DM hindcast: slope = 0.0014 (clears min_slope) but p = 0.42, and the
  # delta-method interval on t* spans centuries and includes the past.
  dm <- c(0.935, 0.902, 0.894, 0.860, 0.921, 0.885, 0.965, 0.901,
          0.938, 0.975, 0.949, 1.038, 0.921, 0.854, 0.902, 0.942)
  h <- ErrorTracer:::.compute_horizon(make_sl_df(dm), threshold = 1.0,
                                      min_slope = 1e-4, min_run = 2L,
                                      max_extrapolation_factor = 5)
  expect_identical(h$type, "lower_bound")
  expect_true(is.na(h$value))
  expect_gt(h$slope, 1e-4)            # magnitude gate alone would have passed
  expect_gt(h$slope_p, 0.05)
  expect_match(h$description, "no trend distinguishable", fixed = TRUE)
})

test_that("a genuinely trending series is still projected", {
  set.seed(42)
  ratios <- 0.30 + 0.04 * seq_len(15) + rnorm(15, 0, 0.05)   # clear upward trend
  h <- ErrorTracer:::.compute_horizon(make_sl_df(ratios), threshold = 1.0,
                                      min_slope = 1e-4, min_run = 2L,
                                      max_extrapolation_factor = 10)
  expect_identical(h$type, "projected")
  expect_lt(h$slope_p, 0.05)
  expect_true(is.finite(h$se_t_star))
})

test_that("projection_alpha = 1 restores unconditional extrapolation", {
  dm <- c(0.935, 0.902, 0.894, 0.860, 0.921, 0.885, 0.965, 0.901,
          0.938, 0.975, 0.949, 1.038, 0.921, 0.854, 0.902, 0.942)
  h <- ErrorTracer:::.compute_horizon(make_sl_df(dm), threshold = 1.0,
                                      min_slope = 1e-4, min_run = 2L,
                                      max_extrapolation_factor = 5,
                                      projection_alpha = 1)
  expect_identical(h$type, "projected")
})
