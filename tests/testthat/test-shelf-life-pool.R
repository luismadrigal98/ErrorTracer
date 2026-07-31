# Pooling shelf lives across forecast origins (ErrorTracer 1.3.1).
#
# R1 (MEE-26-05-477, L195): "determining the forecast limit from one forecast is
# a really bad idea. How does this package handle averaging over repeated
# forecasts and the need to index both start time and lead time?"
#
# The statistical point these tests pin: with right-censored origins, the naive
# mean/median of the OBSERVED horizons is biased low, because the origins that
# get censored are exactly the ones with long horizons. Kaplan-Meier uses the
# censored origins as partial information instead of discarding them.

test_that("censoring is not silently dropped", {
  h <- data.frame(lead     = c(2, 3, 2, 3, 5),
                  censored = c(FALSE, FALSE, FALSE, FALSE, TRUE))
  out <- et_shelf_life_pool(h)
  expect_equal(out$n_origins, 5)
  expect_equal(out$n_observed, 4)
  expect_equal(out$n_censored, 1)
})

test_that("all-censored is reported as a result, not a failure", {
  h <- data.frame(lead = c(16, 13, 11, 8, 5), censored = rep(TRUE, 5))
  out <- et_shelf_life_pool(h)
  expect_equal(out$n_censored, 5)
  expect_true(is.na(out$median_lead))
  expect_match(out$method, "all origins censored", fixed = TRUE)
})

test_that("Kaplan-Meier does not fall below the naive observed median", {
  skip_if_not_installed("survival")
  # Long horizons censored: the naive median of observed values is 2, but the
  # censored origins say the true horizon is at least 9, 10, 11.
  h <- data.frame(lead     = c(2, 2, 3, 9, 10, 11),
                  censored = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE))
  out <- et_shelf_life_pool(h)
  expect_identical(out$method, "kaplan-meier")
  naive <- stats::median(h$lead[!h$censored])
  expect_gte(out$median_lead, naive)
})

test_that("with no censoring the pooled median equals the plain median", {
  skip_if_not_installed("survival")
  h <- data.frame(lead = c(2, 3, 4, 5, 6), censored = rep(FALSE, 5))
  out <- et_shelf_life_pool(h)
  expect_equal(out$median_lead, stats::median(h$lead))
})

test_that("groups are summarised separately", {
  h <- data.frame(
    group    = rep(c("DM", "PP"), each = 5),
    lead     = c(16, 13, 11, 8, 5,  2, 3, 2, 3, 3),
    censored = c(rep(TRUE, 5),      rep(FALSE, 5))
  )
  out <- et_shelf_life_pool(h)
  expect_equal(nrow(out), 2)
  dm <- out[out$group == "DM", ]
  pp <- out[out$group == "PP", ]
  expect_true(is.na(dm$median_lead))     # never reached a horizon
  expect_equal(pp$median_lead, 3)
})

test_that("a list of shelf_life objects is accepted", {
  # Build two minimal et_shelf_life objects by hand: one crossing, one not.
  mk <- function(ratios) {
    df <- data.frame(obs_id = seq_along(ratios), time = seq_along(ratios),
                     ci_width = ratios, plausible_range = 1, ratio = ratios,
                     informative = ratios < 1)
    hor <- ErrorTracer:::.compute_horizon(df, threshold = 1, min_slope = 1e-4,
                                          min_run = 2L)
    structure(df, class = c("et_shelf_life", "data.frame"), horizon = hor)
  }
  crossing <- mk(c(0.5, 0.6, 1.2, 1.3, 1.4))   # sustained from lead 3
  flat     <- mk(rep(0.5, 5))                  # never crosses
  out <- et_shelf_life_pool(list(a = crossing, b = flat))
  expect_equal(out$n_origins, 2)
  expect_equal(out$n_observed, 1)
  expect_equal(out$n_censored, 1)
})

test_that("empty input errors rather than returning a bare NA", {
  expect_error(et_shelf_life_pool(data.frame(lead = numeric(0),
                                             censored = logical(0))),
               "No forecast origins")
})
