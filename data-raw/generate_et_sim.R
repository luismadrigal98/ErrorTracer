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
# Training period : 1995-2014  (20 years per cluster → 40 rows)
# Forecast period : 2015-2019  (5 years per cluster  → 10 rows)
#
# Predictors (all standardised using training-period mean and SD):
#   Tmean  — mean growing-season temperature (°C)
#   PPT    — total growing-season precipitation (mm)
#   SWE    — peak snow water equivalent (mm)
#
# Climate simulation:
#   Tmean: AR(1), phi=0.40, mild warming trend (+0.015°C/yr)
#   PPT:   AR(1), phi=0.40, stationary
#   SWE:   derived from Tmean and PPT (neg. Tmean, pos. PPT)
#
# True data-generating model:
#   Cluster A:  z_diff = 0.50·Tmean − 0.30·PPT + 0.20·SWE + ε,  ε ~ N(0, 0.40)
#   Cluster B:  z_diff = 0.30·Tmean − 0.20·PPT − 0.10·SWE + ε,  ε ~ N(0, 0.50)
# ============================================================

n_train <- 20L
n_fcast <-  5L
n_total <- n_train + n_fcast   # 25 years

years_all   <- 1995:2019
years_train <- 1995:2014
years_fcast <- 2015:2019
train_idx   <- seq_len(n_train)
fcast_idx   <- (n_train + 1L):n_total

# True regression parameters
true_params <- list(
  A = c(intercept =  0.00, Tmean =  0.50, PPT = -0.30, SWE =  0.20, sigma = 0.40),
  B = c(intercept =  0.10, Tmean =  0.30, PPT = -0.20, SWE = -0.10, sigma = 0.50)
)

# ============================================================
# 1. Simulate shared climate time series
# ============================================================
# Tmean: AR(1) with mild warming trend.
#   phi = 0.40 (moderate persistence), drift = 0.015°C/yr
#   This keeps forecast values within ~1.5 SDs of training mean.
Tmean_raw        <- numeric(n_total)
Tmean_raw[1L]    <- 14.5
for (t in 2L:n_total) {
  Tmean_raw[t] <- 0.40 * Tmean_raw[t - 1L] +
                  0.60 * 14.5 +
                  0.015 * t +               # mild warming trend
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

# SWE: physically motivated covariance with Tmean (neg) and PPT (pos)
SWE_raw <- pmax(0,
                -0.55 * Tmean_raw +
                 0.40 * PPT_raw +
                 stats::rnorm(n_total, 50, 8))

# Standardise using training-period statistics only
#   (apply the same transformation to forecast years — correct practice)
tmean_mu <- mean(Tmean_raw[train_idx]);  tmean_sd <- stats::sd(Tmean_raw[train_idx])
ppt_mu   <- mean(PPT_raw[train_idx]);   ppt_sd   <- stats::sd(PPT_raw[train_idx])
swe_mu   <- mean(SWE_raw[train_idx]);   swe_sd   <- stats::sd(SWE_raw[train_idx])

Tmean_std <- (Tmean_raw - tmean_mu) / tmean_sd
PPT_std   <- (PPT_raw   - ppt_mu)   / ppt_sd
SWE_std   <- (SWE_raw   - swe_mu)   / swe_sd

# ============================================================
# 2. Simulate cluster-specific responses
# ============================================================
.make_cluster_data <- function(cluster_name, params, seed_offset) {
  set.seed(111L + seed_offset)

  X      <- cbind(Tmean_std, PPT_std, SWE_std)
  mu_all <- params["intercept"] + X %*% params[c("Tmean", "PPT", "SWE")]
  z_all  <- as.numeric(mu_all) + stats::rnorm(n_total, 0, params["sigma"])

  train_df <- data.frame(
    year       = years_train,
    cluster_id = cluster_name,
    Tmean      = Tmean_std[train_idx],
    PPT        = PPT_std[train_idx],
    SWE        = SWE_std[train_idx],
    z_diff     = z_all[train_idx],
    stringsAsFactors = FALSE
  )

  forecast_df <- data.frame(
    year       = years_fcast,
    cluster_id = cluster_name,
    Tmean      = Tmean_std[fcast_idx],
    PPT        = PPT_std[fcast_idx],
    SWE        = SWE_std[fcast_idx],
    stringsAsFactors = FALSE
  )

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
  # ── Training data: 40 rows (2 clusters × 20 years) ──────────────────────
  train = .sort_df(rbind(cl_A$train, cl_B$train)),

  # ── Forecast predictors: 10 rows (2 clusters × 5 years) ─────────────────
  forecast = .sort_df(rbind(cl_A$forecast, cl_B$forecast)),

  # ── True responses for forecast period (for calibration checks) ──────────
  validation = .sort_df(rbind(cl_A$validation, cl_B$validation)),

  # ── True data-generating parameters (for parameter-recovery checks) ──────
  true_params = true_params,

  # ── Suggested environmental noise SDs for et_predict() ───────────────────
  env_noise = list(Tmean = 0.30, PPT = 0.20, SWE = 0.15),

  # ── Training-period standardisation constants ─────────────────────────────
  standardization = list(
    Tmean = c(mean = tmean_mu, sd = tmean_sd),
    PPT   = c(mean = ppt_mu,   sd = ppt_sd),
    SWE   = c(mean = swe_mu,   sd = swe_sd)
  ),

  description = paste(
    "Simulated allele-frequency-change (z_diff) for two SNP clusters (A, B)",
    "at a mountain plant site. Training: 1995-2014 (20 obs/cluster).",
    "Forecast: 2015-2019 (5 obs/cluster). Three standardised climate",
    "predictors (Tmean, PPT, SWE). True coefficients in et_sim$true_params.",
    "Generated with set.seed(111). See data-raw/generate_et_sim.R."
  )
)

# ============================================================
# 4. Quick sanity checks
# ============================================================
stopifnot(nrow(et_sim$train) == 40L)
stopifnot(nrow(et_sim$forecast) == 10L)
stopifnot(nrow(et_sim$validation) == 10L)
stopifnot(all(c("year","cluster_id","Tmean","PPT","SWE","z_diff") %in% names(et_sim$train)))
stopifnot(!("z_diff" %in% names(et_sim$forecast)))  # predictors only

# Report standardised ranges for quality control
fcast_A <- et_sim$forecast[et_sim$forecast$cluster_id == "A", ]
cat("Forecast Tmean range (std):", round(range(fcast_A$Tmean), 2), "\n")
cat("Forecast PPT range   (std):", round(range(fcast_A$PPT), 2), "\n")
cat("Forecast SWE range   (std):", round(range(fcast_A$SWE), 2), "\n")

# ============================================================
# 5. Save
# ============================================================
dir.create("data", showWarnings = FALSE)
save(et_sim, file = "data/et_sim.rda", compress = "bzip2")
cat("Saved data/et_sim.rda\n")
cat("  Training rows   :", nrow(et_sim$train), "\n")
cat("  Forecast rows   :", nrow(et_sim$forecast), "\n")
cat("  Validation rows :", nrow(et_sim$validation), "\n")
cat("Cluster A true params:", et_sim$true_params$A, "\n")
cat("Cluster B true params:", et_sim$true_params$B, "\n")
