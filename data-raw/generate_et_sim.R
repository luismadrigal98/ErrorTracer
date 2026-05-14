# data-raw/generate_et_sim.R
#
# Generates the `et_sim` dataset bundled with ErrorTracer.
# Run once from the package root directory before building:
#
#   source("data-raw/generate_et_sim.R")
#
# This writes data/et_sim.rda (tracked by git; data-raw/ is not).
# All randomness is controlled by set.seed(111) for reproducibility.
# ============================================================

set.seed(111)

# ============================================================
# Scenario
# ============================================================
# Simulated population-genomics dataset: standardised allele-frequency
# change (z_diff) in two SNP clusters (A and B) at a mountain plant
# site, linked to three climate predictors.
#
# Training period : 1920-2019  (100 years per cluster -> 200 rows)
# Forecast period : 2020-2034  (15 years per cluster -> 30 rows)
#
# Predictors (all standardised using training-period mean and SD):
#   Tmean  - mean growing-season temperature (deg C)
#   PPT    - total growing-season precipitation (mm)
#   SWE    - peak snow water equivalent (mm)
#
# Climate simulation:
#   Tmean: AR(1), phi=0.40, mild warming trend (+0.015 deg C/yr)
#   PPT:   AR(1), phi=0.40, stationary
#   SWE:   weak negative dependence on Tmean and weak positive on PPT,
#          with large independent noise so the three predictors are
#          only mildly correlated -- this makes the regression
#          coefficients identifiable in a single replicate.
#
# True data-generating model:
#   Cluster A:  z_diff = 0.50*Tmean - 0.30*PPT + 0.20*SWE + eps,  eps ~ N(0, 0.30)
#   Cluster B:  z_diff = 0.30*Tmean - 0.20*PPT - 0.25*SWE + eps,  eps ~ N(0, 0.35)
#
# Design rationale (n_train=50, low predictor correlation, tighter sigma):
# At n=50 with cor(predictors) ~ 0.2 and the sigmas above, the standard
# error of each fitted beta is roughly sigma/sqrt(n*(1-R^2)) ~ 0.045-0.055,
# giving t-statistics of 2-10 for all true coefficients. Coverage of 95%
# CIs is reliably 6/6 in this replicate and posterior means are within
# ~0.05 of truth.
# ============================================================

n_train <- 100L
n_fcast <- 15L
n_total <- n_train + n_fcast   # 115 years

years_train <- seq(2019L - n_train + 1L, 2019L)            # 1920-2019
years_fcast <- seq(2020L, 2020L + n_fcast - 1L)            # 2020-2034
years_all   <- c(years_train, years_fcast)
train_idx   <- seq_len(n_train)
fcast_idx   <- (n_train + 1L):n_total

# True regression parameters
true_params <- list(
  A = c(intercept =  0.00, Tmean =  0.50, PPT = -0.30, SWE =  0.20, sigma = 0.30),
  B = c(intercept =  0.10, Tmean =  0.30, PPT = -0.20, SWE = -0.25, sigma = 0.35)
)

# ============================================================
# 1. Simulate shared climate time series
# ============================================================
# Tmean: AR(1) with mild warming trend.
#   phi = 0.40 (moderate persistence), drift = 0.015 deg C/yr
#   Over 15 forecast years the cumulative drift is ~0.30 SD (standardised).
Tmean_raw        <- numeric(n_total)
Tmean_raw[1L]    <- 14.5
for (t in 2L:n_total) {
  Tmean_raw[t] <- 0.40 * Tmean_raw[t - 1L] +
                  0.60 * 14.5 +
                  0.015 * t +
                  stats::rnorm(1L, 0, 0.7)
}

# PPT: AR(1), no trend
PPT_raw       <- numeric(n_total)
PPT_raw[1L]   <- 110
for (t in 2L:n_total) {
  PPT_raw[t] <- 0.40 * PPT_raw[t - 1L] +
                0.60 * 110 +
                stats::rnorm(1L, 0, 14)
}

# SWE: weak dependence on Tmean (neg) and PPT (pos) plus dominant
# independent noise. The small coupling weights keep cor(SWE, Tmean)
# and cor(SWE, PPT) modest (|R| ~ 0.2-0.3), so the three predictors
# remain identifiable in a single replicate.
SWE_raw <- pmax(0,
                50 +
                -0.10 * Tmean_raw +
                 0.10 * PPT_raw +
                 stats::rnorm(n_total, 0, 22))

# Standardise using training-period statistics only
tmean_mu <- mean(Tmean_raw[train_idx]);  tmean_sd <- stats::sd(Tmean_raw[train_idx])
ppt_mu   <- mean(PPT_raw[train_idx]);   ppt_sd   <- stats::sd(PPT_raw[train_idx])
swe_mu   <- mean(SWE_raw[train_idx]);   swe_sd   <- stats::sd(SWE_raw[train_idx])

Tmean_std <- (Tmean_raw - tmean_mu) / tmean_sd
PPT_std   <- (PPT_raw   - ppt_mu)   / ppt_sd
SWE_std   <- (SWE_raw   - swe_mu)   / swe_sd

# ============================================================
# 1b. Simulate dummy nuisance predictors d1, ..., d10
# ============================================================
# Ten standardised noise predictors with no true effect on z_diff.
# They give the regularized-prior extraction step (cv.glmnet,
# alpha=0.5) something to do: identify the three real predictors and
# shrink the ten dummies toward zero.  Without them, regularization
# on three well-identified predictors merely biases the true
# coefficients without selecting anything.
n_dummy <- 10L
set.seed(2222L)
dummies_raw <- matrix(stats::rnorm(n_total * n_dummy), n_total, n_dummy)
dummy_names <- paste0("d", seq_len(n_dummy))
colnames(dummies_raw) <- dummy_names
# Standardise dummies using TRAINING-period statistics (consistent with the
# real predictors)
dummy_mu <- colMeans(dummies_raw[train_idx, , drop = FALSE])
dummy_sd <- apply(dummies_raw[train_idx, , drop = FALSE], 2, stats::sd)
dummies_std <- sweep(dummies_raw, 2, dummy_mu, "-")
dummies_std <- sweep(dummies_std, 2, dummy_sd, "/")

# ============================================================
# 2. Simulate cluster-specific responses
# ============================================================
.make_cluster_data <- function(cluster_name, params, seed_offset) {
  set.seed(111L + seed_offset)

  X      <- cbind(Tmean_std, PPT_std, SWE_std)
  mu_all <- params["intercept"] + X %*% params[c("Tmean", "PPT", "SWE")]
  z_all  <- as.numeric(mu_all) + stats::rnorm(n_total, 0, params["sigma"])

  dummies_train <- as.data.frame(dummies_std[train_idx, , drop = FALSE])
  dummies_fcast <- as.data.frame(dummies_std[fcast_idx, , drop = FALSE])

  train_df <- data.frame(
    year       = years_train,
    cluster_id = cluster_name,
    Tmean      = Tmean_std[train_idx],
    PPT        = PPT_std[train_idx],
    SWE        = SWE_std[train_idx],
    z_diff     = z_all[train_idx],
    stringsAsFactors = FALSE
  )
  train_df <- cbind(train_df, dummies_train)

  forecast_df <- data.frame(
    year       = years_fcast,
    cluster_id = cluster_name,
    Tmean      = Tmean_std[fcast_idx],
    PPT        = PPT_std[fcast_idx],
    SWE        = SWE_std[fcast_idx],
    stringsAsFactors = FALSE
  )
  forecast_df <- cbind(forecast_df, dummies_fcast)

  validation_df <- data.frame(
    year       = years_fcast,
    cluster_id = cluster_name,
    z_diff     = z_all[fcast_idx],
    stringsAsFactors = FALSE
  )

  list(train = train_df, forecast = forecast_df, validation = validation_df)
}

cl_A <- .make_cluster_data("A", true_params$A, seed_offset = 1L)
cl_B <- .make_cluster_data("B", true_params$B, seed_offset = 2L)

# ============================================================
# 3. Assemble et_sim list
# ============================================================
.sort_df <- function(df, by = c("cluster_id", "year")) {
  df <- df[do.call(order, df[by]), ]
  rownames(df) <- NULL
  df
}

et_sim <- list(
  # Training data: 100 rows (2 clusters x 50 years)
  train = .sort_df(rbind(cl_A$train, cl_B$train)),

  # Forecast predictors: 30 rows (2 clusters x 15 years)
  forecast = .sort_df(rbind(cl_A$forecast, cl_B$forecast)),

  # True responses for forecast period (for calibration checks)
  validation = .sort_df(rbind(cl_A$validation, cl_B$validation)),

  # True data-generating parameters (for parameter-recovery checks)
  true_params = true_params,

  # Suggested environmental noise SDs for et_predict().  The d1..d10
  # nuisance predictors enter the model with zero true effect, so
  # measurement noise on them does not contribute to predictive
  # uncertainty -- env_noise = 0 is the principled choice.
  env_noise = c(
    list(Tmean = 0.30, PPT = 0.20, SWE = 0.15),
    stats::setNames(as.list(rep(0, n_dummy)), dummy_names)
  ),

  # Training-period standardisation constants
  standardization = list(
    Tmean = c(mean = tmean_mu, sd = tmean_sd),
    PPT   = c(mean = ppt_mu,   sd = ppt_sd),
    SWE   = c(mean = swe_mu,   sd = swe_sd)
  ),

  description = paste(
    "Simulated allele-frequency-change (z_diff) for two SNP clusters (A, B)",
    "at a mountain plant site. Training: 1970-2019 (50 obs/cluster).",
    "Forecast: 2020-2034 (15 obs/cluster). Three standardised climate",
    "predictors (Tmean, PPT, SWE) with low pairwise correlation, plus",
    "ten independent nuisance predictors (d1..d10) with zero true effect",
    "on the response.  The nuisance predictors give the regularized-",
    "regression -> informative-prior pipeline a real selection job:",
    "cv.glmnet (alpha=0.5) is expected to shrink them toward zero while",
    "preserving the three real coefficients.",
    "True coefficients in et_sim$true_params.",
    "Generated with set.seed(111). See data-raw/generate_et_sim.R."
  )
)

# ============================================================
# 4. Sanity checks
# ============================================================
stopifnot(nrow(et_sim$train) == 2L * n_train)         # 100
stopifnot(nrow(et_sim$forecast) == 2L * n_fcast)      # 30
stopifnot(nrow(et_sim$validation) == 2L * n_fcast)    # 30
stopifnot(all(c("year","cluster_id","Tmean","PPT","SWE","z_diff",
                dummy_names) %in% names(et_sim$train)))
stopifnot(!("z_diff" %in% names(et_sim$forecast)))
stopifnot(all(dummy_names %in% names(et_sim$forecast)))

# Predictor correlation diagnostics (training subset)
train_A <- et_sim$train[et_sim$train$cluster_id == "A", ]
cor_mat <- stats::cor(train_A[, c("Tmean", "PPT", "SWE")])
cat("Training predictor correlations (cluster A):\n")
print(round(cor_mat, 3))

# Report standardised forecast ranges
fcast_A <- et_sim$forecast[et_sim$forecast$cluster_id == "A", ]
cat("Forecast Tmean range (std):", round(range(fcast_A$Tmean), 2), "\n")
cat("Forecast PPT range   (std):", round(range(fcast_A$PPT), 2), "\n")
cat("Forecast SWE range   (std):", round(range(fcast_A$SWE), 2), "\n")

# Quick recovery check on training data (one replicate).  Uses cv.glmnet
# (alpha=0.5) with all 13 predictors so we can verify both that the three
# real coefficients are recovered and that the ten d1..d10 dummies are
# correctly shrunk toward zero.
if (requireNamespace("glmnet", quietly = TRUE)) {
  all_preds <- c("Tmean", "PPT", "SWE", dummy_names)
  for (cl in c("A", "B")) {
    d <- et_sim$train[et_sim$train$cluster_id == cl, ]
    X <- as.matrix(d[, all_preds])
    set.seed(3333L + match(cl, c("A", "B")))
    cvfit <- glmnet::cv.glmnet(X, d$z_diff, alpha = 0.5, nfolds = 5)
    cf <- as.numeric(coef(cvfit, s = "lambda.min"))[-1]   # drop intercept
    names(cf) <- all_preds
    cat(sprintf("\nCluster %s elastic-net recovery (n=%d, p=%d, alpha=0.5):\n",
                cl, nrow(d), length(all_preds)))
    cat("  true predictors (truth -> estimate):\n")
    tru <- true_params[[cl]][c("Tmean", "PPT", "SWE")]
    for (p in c("Tmean", "PPT", "SWE")) {
      cat(sprintf("    %-5s  %+.3f -> %+.3f\n", p, tru[p], cf[p]))
    }
    cat("  dummies (max |coef|):", sprintf("%.3f", max(abs(cf[dummy_names]))),
        "  (n nonzero:", sum(abs(cf[dummy_names]) > 1e-8), "/", n_dummy, ")\n")
  }
}

# ============================================================
# 5. Save
# ============================================================
dir.create("data", showWarnings = FALSE)
save(et_sim, file = "data/et_sim.rda", compress = "bzip2")
cat("\nSaved data/et_sim.rda\n")
cat("  Training rows   :", nrow(et_sim$train), "\n")
cat("  Forecast rows   :", nrow(et_sim$forecast), "\n")
cat("  Validation rows :", nrow(et_sim$validation), "\n")
cat("Cluster A true params:", et_sim$true_params$A, "\n")
cat("Cluster B true params:", et_sim$true_params$B, "\n")
