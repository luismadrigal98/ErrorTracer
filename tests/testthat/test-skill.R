# tests/testthat/test-skill.R
#
# Unit tests for the null-model skill (T6): CRPS helpers, et_skill_score(), and
# the shelf_life() skill gate + default response_scale. No Stan fit needed.

# ── CRPS helpers ─────────────────────────────────────────────────────────────

test_that(".crps_norm matches the closed form for N(0,1) at y=0", {
  # CRPS(N(0,1), 0) = 2*dnorm(0) - 1/sqrt(pi)
  expect_equal(ErrorTracer:::.crps_norm(0, 0, 1),
               2 * dnorm(0) - 1 / sqrt(pi), tolerance = 1e-10)
})

test_that(".crps_sample reduces to |c - y| for a point-mass sample", {
  expect_equal(ErrorTracer:::.crps_sample(rep(5, 200), 3), 2, tolerance = 1e-10)
})

test_that(".crps_sample of a large normal sample ~ .crps_norm", {
  set.seed(1)
  x <- rnorm(20000)
  expect_equal(ErrorTracer:::.crps_sample(x, 0.3),
               ErrorTracer:::.crps_norm(0.3, 0, 1), tolerance = 0.02)
})

# ── mock et_prediction ───────────────────────────────────────────────────────

.mock_skill_pred <- function(draws, newdata, train_y) {
  model <- structure(
    list(data    = data.frame(y = train_y, x = seq_along(train_y)),
         formula = y ~ x),
    class = "et_model")
  structure(list(predictive_draws = draws, newdata = newdata, model = model),
            class = "et_prediction")
}

# ── et_skill_score ───────────────────────────────────────────────────────────

test_that("et_skill_score: an accurate, sharp model beats a wide climatology", {
  set.seed(2)
  n <- 10L; S <- 1000L
  y_true  <- rnorm(n, 0, 0.2)
  train_y <- rnorm(400, 0, 3)                       # wide, uninformative null
  draws   <- sapply(y_true, function(m) rnorm(S, m, 0.4))  # S x n, accurate+sharp
  nd      <- data.frame(x = seq_len(n), t = seq_len(n))
  pred    <- .mock_skill_pred(draws, nd, train_y)

  sk <- et_skill_score(pred, observed = data.frame(y = y_true), response_col = "y",
                       time_col = "t", null = "climatology")
  expect_true(all(c("crps_model", "crps_null", "skill", "skillful") %in% names(sk)))
  expect_gt(mean(sk$skill), 0)                      # model beats null
  expect_true(all(sk$skillful))
  fl <- attr(sk, "forecast_limit")
  expect_equal(fl$type, "lower_bound")              # never loses to null
})

test_that("et_skill_score flags a precise-but-biased forecast as unskillful", {
  set.seed(3)
  n <- 8L; S <- 800L
  train_y <- rnorm(400, 0, 5)                        # wide null that COVERS the obs
  y_true  <- rep(4, n)                               # within the climatology range
  draws   <- sapply(seq_len(n), function(i) rnorm(S, -4, 0.2))  # confident, biased
  nd      <- data.frame(x = seq_len(n), t = seq_len(n))
  pred    <- .mock_skill_pred(draws, nd, train_y)

  sk <- et_skill_score(pred, observed = data.frame(y = y_true), response_col = "y",
                       time_col = "t", null = "climatology")
  expect_lt(mean(sk$skill), 0)                      # loses to null
  expect_equal(attr(sk, "forecast_limit")$type, "observed")
})

test_that("random_walk null derives rw_sd/start from training when unspecified", {
  set.seed(4)
  n <- 6L; S <- 500L
  train_y <- cumsum(rnorm(50))                       # a random walk itself
  y_true  <- rnorm(n, tail(train_y, 1), 1)
  draws   <- sapply(y_true, function(m) rnorm(S, m, 1))
  pred    <- .mock_skill_pred(draws, data.frame(x = seq_len(n)), train_y)
  sk <- et_skill_score(pred, observed = data.frame(y = y_true),
                       response_col = "y", null = "random_walk")
  expect_equal(nrow(sk), n)
  expect_true(all(is.finite(sk$crps_null)))
})

# ── sustained forecast limit (regression: isolated dip must not set it) ──────
#
# Same defect class as shelf_life()'s first-exceedance bug: skill hovering near
# zero crosses by chance, and a bare first-crossing rule promotes one lead to
# "the forecast limit". Found on the real Kyoto/climatology comparison, whose
# mean skill is comfortably positive across the window.

.mock_skill_pred_with_dip <- function(bad_at, n = 10L, S = 1000L, seed = 11L) {
  set.seed(seed)
  train_y <- rnorm(400, 0, 3)                    # wide climatology
  y_true  <- rnorm(n, 0, 0.2)
  draws   <- sapply(seq_len(n), function(i) {
    if (i %in% bad_at) rnorm(S, y_true[i] + 12, 0.2)  # confident and wrong
    else               rnorm(S, y_true[i], 0.4)       # accurate and sharp
  })
  list(pred = .mock_skill_pred(draws, data.frame(x = seq_len(n), t = seq_len(n)),
                               train_y),
       y_true = y_true)
}

test_that("an isolated non-skillful lead does not set the forecast limit", {
  f  <- .mock_skill_pred_with_dip(bad_at = 4L)
  sk <- et_skill_score(f$pred, observed = data.frame(y = f$y_true),
                       response_col = "y", time_col = "t", null = "climatology")
  expect_false(sk$skillful[4])
  expect_true(all(sk$skillful[-4]))

  fl <- attr(sk, "forecast_limit")
  expect_equal(fl$type, "lower_bound")   # NOT "observed" at t = 4
  expect_true(is.na(fl$value))
  expect_equal(fl$n_below, 1L)
  expect_equal(fl$first_below, 4)        # the dip is reported, not hidden
  expect_equal(fl$min_run, 2L)
})

test_that("min_run = 1 restores the first-crossing behaviour", {
  f  <- .mock_skill_pred_with_dip(bad_at = 4L)
  sk <- et_skill_score(f$pred, observed = data.frame(y = f$y_true),
                       response_col = "y", time_col = "t", null = "climatology",
                       min_run = 1L)
  fl <- attr(sk, "forecast_limit")
  expect_equal(fl$type, "observed")
  expect_equal(fl$value, 4)
})

test_that("a sustained run of non-skillful leads does set the limit", {
  f  <- .mock_skill_pred_with_dip(bad_at = c(6L, 7L, 8L))
  sk <- et_skill_score(f$pred, observed = data.frame(y = f$y_true),
                       response_col = "y", time_col = "t", null = "climatology")
  fl <- attr(sk, "forecast_limit")
  expect_equal(fl$type, "observed")
  expect_equal(fl$value, 6)              # first lead of the run, not the last
  expect_equal(fl$n_below, 3L)
})

test_that("min_run is honoured for grouped predictions", {
  f1 <- .mock_skill_pred_with_dip(bad_at = 4L, seed = 11L)
  f2 <- .mock_skill_pred_with_dip(bad_at = c(6L, 7L), seed = 12L)
  plist <- structure(
    list(predictions = list(A = f1$pred, B = f2$pred), grouping = "sp"),
    class = "et_prediction_list")
  obs <- rbind(data.frame(sp = "A", y = f1$y_true),
               data.frame(sp = "B", y = f2$y_true))
  sk <- et_skill_score(plist, observed = obs, response_col = "y",
                       time_col = "t", null = "climatology")
  lims <- attr(sk, "forecast_limit_by_group")
  expect_equal(lims$A$type, "lower_bound")   # isolated dip
  expect_equal(lims$B$type, "observed")      # sustained run
  expect_equal(lims$B$value, 6)
})

# ── shelf_life() skill gate + default response_scale ─────────────────────────

.mock_shelf_pred <- function(widths, train_y, times = seq_along(widths)) {
  n <- length(widths)
  ci <- data.frame(row_id = seq_len(n), ci_level = 0.90,
                   lower = -widths / 2, median = 0, upper = widths / 2,
                   width = widths, stringsAsFactors = FALSE)
  model <- structure(list(data = data.frame(y = train_y, x = seq_along(train_y)),
                          formula = y ~ x), class = "et_model")
  structure(list(credible_intervals = ci,
                 newdata = data.frame(x = seq_len(n), t = times),
                 model = model),
            class = "et_prediction")
}

test_that("shelf_life skill gate turns precise-but-unskillful periods uninformative", {
  pred <- .mock_shelf_pred(widths = rep(0.2, 6), train_y = rnorm(100, 0, 2))
  # all precise (narrow CI vs range ~ several units)
  sl_prec <- shelf_life(pred, response_scale = c(-3, 3), ci_level = 0.90)
  expect_true(all(sl_prec$informative))

  skill_df <- data.frame(obs_id = 1:6, skillful = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
  sl_gate  <- shelf_life(pred, response_scale = c(-3, 3), ci_level = 0.90,
                         skill = skill_df)
  expect_equal(sl_gate$informative, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
})

test_that("shelf_life defaults response_scale to the training response range", {
  train_y <- c(-4, 10, rnorm(50))
  pred <- .mock_shelf_pred(widths = rep(0.5, 4), train_y = train_y)
  sl <- suppressMessages(shelf_life(pred, ci_level = 0.90))   # no response_scale
  expect_equal(unique(sl$plausible_range), diff(range(train_y)), tolerance = 1e-8)
})
