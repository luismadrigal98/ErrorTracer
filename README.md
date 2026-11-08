# ErrorTracer: Bayesian Error Propagation and Forecast Uncertainty Decomposition

## Overview

<img src="images/logo.png" align="right" height="200" alt="ErrorTracer logo" />

**ErrorTracer** provides a complete pipeline for ecological and genomic forecasting with climate or environmental covariates. It bridges the gap between regularized regression and fully Bayesian workflows in three steps:

1. **Extract informed priors** from a fitted regularized or standard model (`glmnet`, `lm`, `glm`, `ranger`) to set data-informed prior scales (default `shrinkage = "zero"`, centered at zero to avoid prior/likelihood double-use, or `shrinkage = "estimate"` to carry point estimates forward as prior means; `et_priors_split()` splits the sample so selection and likelihood share no rows) — or bring an existing `brms` fit straight into the pipeline via `as_et_model()`.
2. **Decompose forecast uncertainty** into independent sources — parameter uncertainty, environmental/covariate measurement or trajectory uncertainty, and residual process variance — plus a group-level component when predicting for new random-effect levels and a temporal-autocorrelation component whenever the model formula carries an `ar()`, `ma()`, `arma()`, `cosy()`, `unstr()`, `sar()`, or `car()` term (alongside order-independent global Sobol indices via `et_sobol()`).
3. **Quantify forecast shelf life** — the time horizon at which a forecast's credible interval becomes wider than the plausible response scale, rendering it uninformative, with optional CRPS skill gating (`et_skill_score()`) and Kaplan–Meier pooling across rolling forecast origins (`et_shelf_life_pool()`).

The package is designed for ecological and genomic time-series forecasting in which regression models are fitted or refitted in `brms` and predictions must be accompanied by a principled uncertainty budget. The API is fully general — response, predictors, and grouping are formula-driven.

---

## Key Features

| Feature | Description |
|---|---|
| **Regularized → Bayesian prior pipeline** | `extract_priors()` converts `glmnet`, `lm`, `glm`, or `ranger` fits into `brms`-compatible priors (`et_priors_split()` for disjoint data-splitting); or pass an existing `brmsfit` directly via `as_et_model()` |
| **Uncertainty decomposition** | Partitions measured posterior predictive variance (`total_var`) into parameter, environmental, and residual components (plus `group_var` for new random-effect levels and signed `temporal_var` for autocorrelation models), with global Sobol indices via `et_sobol()` |
| **Forecast shelf life & pooling** | `shelf_life()` quantifies when a forecast becomes uninformative (with persistence and trend-significance gates); `et_shelf_life_pool()` pools horizons across origins via Kaplan–Meier |
| **Skill scoring & calibration** | `et_skill_score()` scores forecasts against climatology or random-walk nulls (CRPS); `et_calibrate()` and `et_pit()` assess interval coverage and PIT uniformity |
| **Group-level fitting** | One function call fits and predicts across all groups (e.g., SNP clusters, species, sites) |
| **Rich visualizations** | Eight manuscript-ready `ggplot2`-based plotting functions (`et_plot_*()`), all returning customizable objects |
| **Full S3 class system** | Consistent `print`, `summary`, and `plot`-family methods across all objects |

---

## Installation

Install the development version from GitHub using `remotes`:

```r
# install.packages("remotes")
remotes::install_github("madrigalrocalj/ErrorTracer")
```

ErrorTracer requires a Stan backend. Install either `cmdstanr` (recommended) or `rstan` separately:

```r
# Option A: cmdstanr (recommended)
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

# Option B: rstan
install.packages("rstan")
```

---

## Core Workflow

### Step 1 — Extract informed priors from a regularized model

`extract_priors()` is an S3 generic that dispatches on the class of the fitted model. With the default `shrinkage = "zero"`, coefficient priors are centered at zero with SDs scaled by coefficient magnitude (`multiplier * |coef|`, floored at `min_sd`) — a regularizing, ridge-style prior that avoids reusing point estimates in the likelihood. Set `shrinkage = "estimate"` to use the fitted coefficients as informative prior means.

```r
library(ErrorTracer)
library(glmnet)

# --- Elastic net (cv.glmnet) ---
cv_fit <- cv.glmnet(x_train, y_train, alpha = 0.5)

prior_spec <- extract_priors(
  model              = cv_fit,
  lambda             = "lambda.min",  # or "lambda.1se"
  multiplier         = 2.0,           # prior SD = multiplier * |coef|
  min_sd             = 0.1,           # floor for prior SD
  shrinkage          = "zero",        # default: mean 0; or "estimate" for prior mean = coef
  intercept_prior_sd = NULL,          # NULL lets brms use its data-aware default
  sigma_prior_scale  = 1.0            # half-Cauchy scale for sigma
)

print(prior_spec)
# ErrorTracer prior specification
#   Method      : glmnet
#   Predictors  : 5
#   Multiplier  : 2
#   Min SD      : 0.1
#   Shrinkage   : zero (prior means at 0; regularizing)
#   Coefficients:
#     Tmean                 prior mean =   0.0000  sd =   0.6482  (coef =   0.3241)
#     PPT                   prior mean =   0.0000  sd =   0.2174  (coef =  -0.1087)
#     ...
```

**Supported model classes and their prior construction:**

| Class | Prior mean (`shrinkage = "zero"`, default) | Prior mean (`shrinkage = "estimate"`) | Prior SD |
|---|---|---|---|
| `cv.glmnet` / `glmnet` | `0` (for non-zero coefs) | Non-zero coefficients | `multiplier × |coef|` (floor: `min_sd`) |
| `lm` | `0` | Coefficient estimates | `multiplier × |coef|` (`"zero"`) or `multiplier × SE` (`"estimate"`, floor: `min_sd`) |
| `glm` | `0` | Coefficient estimates | `multiplier × |coef|` (`"zero"`) or `multiplier × SE` (`"estimate"`, floor: `min_sd`) |
| `ranger` (experimental) | `0` (direction unknown) | `0` (direction unknown) | `multiplier × normalised_importance` (floor: `min_sd`) |

For `ranger`, only variables with positive permutation importance are included. Because random forests do not produce signed coefficients or GLM-scale units, `extract_priors.ranger()` is experimental and centers priors at zero scaled by normalized importance.

> **Note on intercept handling:** The first-stage intercept is excluded from the coefficient priors, and `intercept_prior_sd = NULL` (default) lets `brms` use its data-aware default (`student_t(3, median(y), mad(y))`). Supplying a numeric `intercept_prior_sd` sets `normal(0, intercept_prior_sd)`, which is only appropriate when the response is centered.

---

### Step 2 — Fit the Bayesian model

`et_fit()` wraps `brms::brm()` and attaches the prior specification and training data for downstream use.

```r
fit <- et_fit(
  formula       = z_diff ~ Tmean + PPT + SWE,
  data          = train_df,
  priors        = prior_spec,       # from extract_priors(); or NULL for brms defaults
  chains        = 4L,
  iter          = 2000L,
  warmup        = 1000L,
  cores         = 4L,
  seed          = 42L,
  adapt_delta   = 0.95,
  max_treedepth = 12L
)

print(fit)
# ErrorTracer model (et_model)
#   Formula : z_diff ~ Tmean + PPT + SWE
#   n obs   : 48
#   Chains  : 4  Iter: 2000  Warmup: 1000
#   Priors  : informed (glmnet, 3 predictors)
#   Rhat max: 1.002

summary(fit)
```

#### Group-level fitting

Pass a `grouping` column name to fit one model per group and receive an `et_model_list`:

```r
# Fits one brms model per SNP cluster
fit_list <- et_fit(
  formula  = z_diff ~ Tmean + PPT + SWE,
  data     = all_clusters_df,
  priors   = prior_spec,        # single spec applied to all groups
  grouping = "cluster_id",      # or pass a named list of et_prior_spec for per-group priors
  chains   = 4L,
  iter     = 2000L
)

print(fit_list)
# ErrorTracer grouped model list (et_model_list)
#   Grouping : cluster_id
#   Formula  : z_diff ~ Tmean + PPT + SWE
#   Groups   : 12
#   Fitted   : 12 / 12
```

#### Bringing your own fit

Steps 1 and 2 are optional. A model you fitted yourself with `brms` runs
through the rest of the pipeline unchanged — `et_predict()` and `et_diagnose()`
accept a `brmsfit` directly, and `as_et_model()` gives you the `et_model`
object if you want it:

```r
my_fit <- brms::brm(z_diff ~ Tmean + PPT, data = train_df)

pred <- et_predict(my_fit, newdata = forecast_df,
                   env_noise = list(Tmean = 0.5, PPT = 0.1))

# Or wrap it explicitly, e.g. to attach a prior spec for the
# prior-vs-posterior plots, or to supply training data the formula did not use:
model <- as_et_model(my_fit, prior_spec = prior_spec)
```

The formula, family, training data and sampling seed are read off the fit, so a
prediction from a wrapped model is as reproducible as one from `et_fit()`.

Posterior draws on their own are not enough: a bare `stanfit` or `CmdStanMCMC`
records nothing about which parameters are coefficients, what the link function
is, or which rows a forecast continues. To bring draws from a hand-written Stan
program in, give `brms` the formula, family and data once with
`brm(..., empty = TRUE)` and attach the draws to the result — the error message
spells out the four lines.

---

### Step 3 — Predict with uncertainty propagation

`et_predict()` draws from the posterior predictive distribution, propagates environmental measurement noise (or pathwise driver ensembles via `env_ensemble`) through the model, and computes credible intervals. By default (`seed = NULL`), predictions are seeded from the fitted model's own seed so results are reproducible from the model alone.

```r
pred <- et_predict(
  model     = fit,
  newdata   = forecast_df,
  env_noise = list(Tmean = 0.5, PPT = 10.0, SWE = 5.0),
  # env_noise can also be a single fraction of each predictor's SD:
  # env_noise = 0.10   →  10% noise applied uniformly
  n_draws   = 2000L,
  ci_levels = c(0.5, 0.8, 0.90, 0.95)
)

print(pred)
# ErrorTracer prediction (et_prediction)
#   Observations : 24
#   Draws        : 2000
#   CI levels    : 0.5, 0.8, 0.9, 0.95
#   Mean var decomposition (across observations):
#     Parameter    : 0.0412
#     Environmental: 0.0083
#     Residual     : 0.1204
#     Total        : 0.1699
```

The returned `et_prediction` object contains:

| Slot | Description |
|---|---|
| `posterior_predict` | `[n_draws × n_obs]` matrix — unperturbed posterior predictive draws |
| `posterior_linpred` | `[n_draws × n_obs]` matrix — linear predictor draws on the link scale (parameter uncertainty only) |
| `lp_perturbed` | `[n_perturb × n_obs]` matrix — linear predictor on the link scale under environmentally perturbed inputs |
| `predictive_draws` | `[n_draws × n_obs]` matrix — the predictive draws used for `credible_intervals`, `et_pit()`, and `et_skill_score()` (includes environmental uncertainty when `include_env_in_ci = TRUE`) |
| `sigma_draws` | Numeric vector of posterior sigma draws (`NA` for families without sigma) |
| `credible_intervals` | `data.frame` with columns `row_id, ci_level, lower, median, upper, width` |
| `decomposition` | `data.frame` with `obs_id, total_var, param_var, env_var, v_env_mcse, residual_var` (plus `group_var` for new random-effect levels, and `temporal_var` for autocorrelation models or `closure_gap` for iid ones) |
| `newdata` | The input forecast data |
| `model` | Reference to the source `et_model` |
| `env_noise`, `env_cov`, `env_dist`, `env_ensemble` | Environmental perturbation / ensemble specifications used |

---

### Step 4 — Decompose uncertainty

`decompose_uncertainty()` extracts the variance decomposition table from a prediction object:

```r
decomp <- decompose_uncertainty(pred)
head(decomp)
#   obs_id total_var param_var env_var v_env_mcse residual_var closure_gap
# 1      1    0.1665    0.0389  0.0071     0.0012       0.1204      0.0001
# 2      2    0.1741    0.0441  0.0098     0.0015       0.1204     -0.0002
# ...
```

**Decomposition method:**

| Component | Estimator |
|---|---|
| `param_var` | `Var[μ]` across draws on the response scale (population level when predicting for new group levels) — coefficient uncertainty |
| `group_var` | `Var` contributed by the group-level distribution — present **only** when predicting for a *new* level of a random effect (see below) |
| `env_var` | `Var[μ_perturbed] − Var[μ_unperturbed]` — additional response-scale variance from predictor noise |
| `v_env_mcse` | Monte Carlo standard error of `env_var` |
| `residual_var` | `E[V(μ)]` across posterior draws (`E[σ²]` for Gaussian) — biological process noise |
| `temporal_var` | Signed remainder `total_var − (param_var + env_var + residual_var)` — autocorrelation-induced accumulation, present only for models carrying a correlation term |
| `closure_gap` | Signed remainder `total_var − (param_var + env_var + residual_var)` for iid models — Monte Carlo closure check |
| `total_var` | **`Var[pp_env]`**, the measured variance of the environmentally inclusive posterior predictive draws |

Since 1.4.0 `total_var` is **measured**, not defined: it is `Var[posterior_predict]`, the spread of the distribution the reported interval is actually drawn from. For a model carrying a correlation term `temporal_var` is the remainder and the four components account for it exactly; for an iid model the shortfall is returned as `closure_gap` rather than folded into a component, because a gap much larger than the Monte Carlo error of its parts means a mis-estimated component rather than noise. A budget whose total is *defined* as the sum of its parts cannot report when the parts fail to account for it. (Before 1.3.0 the total was computed independently and `env_var` sat outside it; 1.3.x then defined it as the component sum. Decomposition tables from either should be recomputed rather than compared.)

**Autocorrelation.** For models whose formula carries `ar()`, `ma()`, `arma()`, `cosy()`, `unstr()`, `sar()` or `car()`, a `temporal_var` column appears. It captures the autocorrelation-induced spread that `brms::posterior_predict()` accumulates iteratively beyond a single innovation, and `residual_var` is then read as the *innovation* variance rather than the stationary marginal variance. Since 1.4.0 `temporal_var` is **signed**: it is a covariance term, and a negative value is meaningful, since autocorrelated errors can reduce predictive variance relative to an independent-innovation baseline. The `pmax(., 0)` guard remains on `env_var`, where a difference of two estimated variances is noisy around a true zero, and it is instrumented: `v_env_mcse` is reported alongside and the floor warns when it fires beyond twice that error, rather than silently absorbing a systematic negative.

**Hierarchical models.** Predicting for a group the fit has *seen* adds no component — that group's effect is an estimated parameter and its uncertainty is already inside `param_var`. Predicting for a **new** level integrates over the group-level distribution and adds `group_var`, carrying the between-group variance. `et_predict()` detects new levels automatically and reports which budget it computed; the two are genuinely different objects and should not be compared directly.

---

### Step 5 — Forecast shelf life

`shelf_life()` computes the ratio of credible interval width to the plausible response range at each forecast time point, and flags when the forecast becomes uninformative:

```r
# For arcsin-sqrt transformed responses, derive range from training data:
#   response_scale = range(train_df$z_diff)
# For bounded responses, supply explicit bounds:
#   response_scale = c(lower, upper)
sl <- shelf_life(
  predictions              = pred,
  response_scale           = range(train_df$z_diff),  # defaults to the training range
  ci_level                 = 0.90,   # must be present in the et_prediction object
  threshold                = 1.0,    # CI/range ratio above which forecast is uninformative
  time_col                 = "year", # column in newdata to use as time axis
  min_run                  = 2L,     # consecutive periods required before a crossing counts
  projection_alpha         = 0.05,   # a projected horizon needs a significantly positive slope
  max_extrapolation_factor = 10      # cap on linear projection beyond observed window
)

print(sl)
# ErrorTracer shelf life analysis
#   Observations     : 24
#   Plausible range  : 2
#   Threshold        : 1
#   Informative      : 18 / 24
#   Mean CI/range    : 0.847
#   Max CI/range     : 1.231
#   Shelf life       : ~2026 (observed -- threshold first exceeded)
#   Detail           : Threshold (ratio >= 1) sustained (>= 2 consecutive periods) from 2026.

as.data.frame(sl)
#   obs_id  time ci_width plausible_range  ratio informative
# 1      1  2024    0.841               2  0.421        TRUE
# 2      2  2025    1.103               2  0.552        TRUE
# 3      3  2026    2.461               2  1.231       FALSE
```

A `threshold` of 1.0 (default) means the forecast is considered uninformative when the CI spans the entire plausible response range. Lower thresholds (e.g., 0.8) impose stricter requirements.

**Two gates guard against reporting noise as a horizon** (both added in 1.3.1, both changing the default output):

- `min_run` (default `2`) requires the crossing to **persist**. A ratio hovering near the threshold without trending will exceed it once by chance, and the previous first-exceedance rule promoted that single excursion to a horizon — on a real 16-period series the reported horizon moved between four different answers under sampler settings alone.
- `projection_alpha` (default `0.05`) requires a projected horizon to rest on a slope that is *significantly* positive. `min_slope_for_projection` is a magnitude gate only, so a flat, noisy trend could be extrapolated to a confident-looking crossing time whose interval spanned centuries.

Every horizon carries `n_exceedances`, `frac_exceedance`, `first_exceedance` and `min_run` (plus `slope`/`slope_p` in projection mode) so an isolated excursion is visible rather than silent. Set `min_run = 1` and `projection_alpha = 1` to reproduce pre-1.3.1 behaviour. **Any horizon computed with an earlier release should be recomputed, not trusted.**

---

### Step 6 — Calibrate against validation data

`et_calibrate()` computes observed coverage probability at each nominal CI level. A well-calibrated model produces coverage that matches the nominal level:

```r
cal <- et_calibrate(
  predictions  = pred,
  observed     = validation_df,  # same n_rows as newdata; matched positionally
  response_col = "z_diff",       # inferred from formula if omitted
  ci_levels    = c(0.5, 0.8, 0.90, 0.95)
)

print(cal)
#   ci_level nominal observed_coverage n_obs calibration_error
# 1     0.50    0.50             0.542    24             0.042
# 2     0.80    0.80             0.792    24             0.008
# 3     0.90    0.90             0.875    24             0.025
# 4     0.95    0.95             0.958    24             0.008
```

---

### Step 7 — Diagnose convergence

`et_diagnose()` reports Rhat, effective sample size ratios, divergent transitions, and LOO-CV in a single call:

```r
diag <- et_diagnose(fit, loo = TRUE)

# Convergence
diag$convergence$rhat_max     # should be < 1.05
diag$convergence$neff_all_ok  # all Neff ratios > 0.1
diag$convergence$n_divergences

# LOO-CV
diag$loo$elpd_loo
diag$loo$n_bad_pareto_k  # Pareto k > 0.7 indicate influential observations
```

For grouped models, `et_diagnose()` returns a `per_group` list and a `summary` data.frame with one row per group.

---

## Forecast evaluation

Shelf life is a statement about **precision** on the response scale. A forecast can be
precise and still wrong, so the package pairs it with an accuracy criterion and a
full-distribution calibration check.

```r
# Null-relative forecast limit (CRPS). Two nulls, hard in opposite regimes:
#   "climatology"  (default) -- the unconditional TRAINING-response distribution:
#                   what a forecaster with no covariates would have issued. Hard to
#                   beat for a mean-reverting series. Needs no external data.
#   "random_walk"            -- persistence: the last training value carried forward
#                   with SD growing as sd(diff(train)) * sqrt(lead). Hard to beat for
#                   a near-unit-root series.
# Report BOTH: which null is demanding is a property of the series, not the model.
sk_clim <- et_skill_score(pred, observed = valid_df, time_col = "year",
                          null = "climatology")
sk_rw   <- et_skill_score(pred, observed = valid_df, time_col = "year",
                          null = "random_walk")
attr(sk_clim, "forecast_limit")   # $value, $type, $n_below, $first_below, $min_run

# Gate the shelf life on skill, so a precise-but-biased forecast cannot pass:
sl_gated <- shelf_life(pred, response_scale = range(train_df$z_diff),
                       ci_level = 0.90, skill = sk_clim)

# Probability integral transform -- the continuous-response analogue of the rank
# histogram. Uniform = calibrated; U-shape = over-confident; tilt = biased.
pit <- et_pit(pred, observed = valid_df)
et_plot_pit(pit)

# Order-independent (Sobol) decomposition, as a cross-check on the additive budget's
# ordering: reports first-order parameter and environmental shares plus their interaction.
sob <- et_sobol(pred, seed = 1L)
```

`et_skill_score()` takes the same `min_run` as `shelf_life()` (default `2`) for the same
reason: skill hovering near zero crosses by chance, and a bare first-crossing rule turns one
unlucky lead into a forecast limit.

## Pooling across forecast origins

A horizon read off a **single** forecast origin inherits whatever was idiosyncratic about the
final training years. `et_shelf_life_pool()` pools horizons from repeated origins, treating
origins whose data ran out before the forecast degraded as **right-censored** rather than
dropping them — dropping them biases the pooled horizon downward, because those are exactly
the long horizons.

```r
# horizons: a data.frame with one row per forecast origin and columns
#   lead     -- the horizon in lead time
#   censored -- TRUE if the forecast was still informative when the data ran out
#   group    -- optional; pools separately per group (e.g. per species)
pooled <- et_shelf_life_pool(horizons, conf_level = 0.95)   # Kaplan-Meier via {survival}
```

## Avoiding the post-selection problem

When the first-stage model is an elastic net, variable selection and the likelihood normally
see the same rows, so reported coverage is conditional on the selected set being correct.
`et_priors_split()` removes the issue by construction:

```r
sp <- et_priors_split(train_df, prior_fun = my_enet_fun, prop = 0.5, seed = 1L)
fit <- et_fit(z_diff ~ Tmean + PPT + SWE, data = sp$fit_data, priors = sp$priors)
```

The cost is sample size: the likelihood sees fewer rows and the posterior is correspondingly
wider. Split when selection stability is in doubt; don't when the predictor set is fixed a
priori and there is no selection to protect against.

---

## Visualization

All plotting functions return `ggplot2` objects that can be further customized with standard `ggplot2` additions (`+ theme(...)`, `+ labs(...)`, etc.).

They deliberately set no figure title or subtitle. These plots are meant for manuscripts, where the caption is typeset separately and a title baked into the image duplicates it. Add one when you want it:

```r
et_plot_pit(pit) + ggplot2::labs(title = "Posterior-predictive calibration")
```

```r
# Stacked bar: proportional contribution of each variance component
et_plot_decomposition(decomp, proportional = TRUE)

# Line chart: CI width / plausible range over time, with uninformative threshold
et_plot_shelf_life(sl, show_ratio = TRUE)

# Calibration: observed vs. nominal coverage (points on the 1:1 diagonal = well-calibrated)
et_plot_calibration(cal)

# PIT histogram: uniformity check with 95% consistency band (per-panel when grouped)
et_plot_pit(pit)

# Fan chart: nested CI ribbons with optional observed values overlaid
et_plot_forecast(pred, observed = validation_df,
                 response_col = "z_diff", time_col = "year")

# Density overlay: prior vs. posterior for each regression coefficient
et_plot_prior_posterior(fit, max_preds = 8L)

# Forest plot: Bayesian 95% CI (blue) vs. regularized coefficient (red ×)
et_plot_coefficients(fit)

# Sensitivity profile: how env_share, CI/range ratio, or horizon responds to noise scale
sens <- et_sensitivity_profile(fit, newdata = forecast_df,
                               response_scale = range(train_df$z_diff),
                               time_col = "year")
et_plot_sensitivity(sens, show = "env_share")
```

---

## S3 Class System

| Class | Created by | Methods |
|---|---|---|
| `et_prior_spec` | `extract_priors()` | `print` |
| `et_model` | `et_fit()`, `as_et_model()` | `print`, `summary` |
| `et_model_list` | `et_fit(..., grouping = ...)` | `print`, `summary` |
| `et_prediction` | `et_predict()` | `print` |
| `et_prediction_list` | `et_predict()` on `et_model_list` | `print` |
| `et_shelf_life` | `shelf_life()` | `print`, `summary` |
| `et_sensitivity` | `et_sensitivity_profile()` | `print` |
| `et_sobol` | `et_sobol()` | `print` |

---

## Full Example: Allele Frequency Forecasting

```r
library(ErrorTracer)
library(glmnet)

# --- Prepare data ---
# train_df: historical years with columns z_diff, Tmean, PPT, SWE
# forecast_df: future years with predictor values only

# 1. Fit elastic net
x_mat <- as.matrix(train_df[, c("Tmean", "PPT", "SWE")])
cv_fit <- cv.glmnet(x_mat, train_df$z_diff, alpha = 0.5)

# 2. Extract informed priors
prior_spec <- extract_priors(cv_fit, multiplier = 2.0, min_sd = 0.1)

# 3. Fit Bayesian model
fit <- et_fit(
  formula = z_diff ~ Tmean + PPT + SWE,
  data    = train_df,
  priors  = prior_spec,
  chains  = 4L, iter = 2000L, cores = 4L
)

# 4. Predict with environmental noise propagation
pred <- et_predict(
  model     = fit,
  newdata   = forecast_df,
  env_noise = list(Tmean = 0.5, PPT = 10.0, SWE = 5.0),
  ci_levels = c(0.5, 0.80, 0.90, 0.95)
)

# 5. Decompose uncertainty
decomp <- decompose_uncertainty(pred)

# 6. Assess forecast shelf life
# For arcsin-sqrt responses use the observed training range; for bounded
# responses supply explicit bounds.
sl <- shelf_life(
  pred,
  response_scale           = range(train_df$z_diff),
  ci_level                 = 0.90,
  threshold                = 1.0,
  time_col                 = "year",
  max_extrapolation_factor = 10
)

# 7. Calibrate (if validation data available)
cal <- et_calibrate(pred, observed = valid_df, response_col = "z_diff")

# 8. Diagnose
et_diagnose(fit)

# 9. Visualize
library(ggplot2)
et_plot_decomposition(decomp)
et_plot_shelf_life(sl)
et_plot_calibration(cal)
et_plot_forecast(pred, observed = valid_df,
                 response_col = "z_diff", time_col = "year")
et_plot_prior_posterior(fit)
et_plot_coefficients(fit)
```

---

## Group-Level Analysis Example

```r
# Fit one model per SNP cluster, then predict and summarize shelf life by cluster
fit_list <- et_fit(
  formula  = z_diff ~ Tmean + PPT + SWE,
  data     = all_df,
  priors   = prior_spec,
  grouping = "cluster_id",
  chains = 4L, iter = 2000L, cores = 4L
)

pred_list <- et_predict(
  fit_list,
  newdata   = forecast_df,
  env_noise = list(Tmean = 0.5, PPT = 10.0),
  ci_levels = c(0.50, 0.90, 0.95)
)

# Per-group decomposition and shelf life
decomp_all <- decompose_uncertainty(pred_list)
sl_all     <- shelf_life(pred_list,
                         response_scale           = range(train_df$z_diff),
                         ci_level                 = 0.90,
                         time_col                 = "year",
                         max_extrapolation_factor = 10)

# Per-group diagnostics summary
diag_list  <- et_diagnose(fit_list)
diag_list$summary  # data.frame: rhat_ok, neff_ok, n_divergences, elpd_loo per cluster

# Faceted shelf life plot
et_plot_shelf_life(sl_all) +
  ggplot2::facet_wrap(~ group, ncol = 4)
```

---

## Comparison with Existing R Packages

| Package | What it does | What is missing |
|---|---|---|
| `brms` | Bayesian regression | No shelf life, no decomposition, no enet→prior pipeline |
| `propagate` | Analytical error propagation | No posterior sampling, no Bayesian framework |
| `rstanarm` | Bayesian regression | Same gaps as `brms` |
| `forecast` / `fable` | Time series forecasting | No ecological/genomic data structures, no decomposition |
| `INLA` | Spatial/temporal Bayes | No shelf life concept, very different user base |

**ErrorTracer's niche:** the full pipeline from regularized feature selection → informed Bayesian refitting → structured uncertainty decomposition → forecast shelf life, designed for ecology and genomics.

---

## Package Structure

```
ErrorTracer/
├── R/
│   ├── priors.R              # extract_priors() — dispatch for glmnet, lm, glm, ranger
│   ├── priors_split.R        # et_priors_split() — disjoint prior/likelihood sample split
│   ├── fit.R                 # et_fit(), as_et_model() — Bayesian fitting & brmsfit wrapping
│   ├── predict.R             # et_predict() — posterior prediction with noise/ensemble propagation
│   ├── decompose.R           # decompose_uncertainty() — param / group / env / residual / temporal split
│   ├── sobol.R               # et_sobol() — global Sobol variance-based sensitivity indices
│   ├── sensitivity.R         # et_sensitivity_profile() — noise-grid sensitivity profiling
│   ├── shelf_life.R          # shelf_life() — forecast horizon analysis
│   ├── shelf_life_rolling.R  # et_shelf_life_pool() — Kaplan-Meier pooling across forecast origins
│   ├── skill.R               # et_skill_score() — null-model CRPS forecast skill & forecast limit
│   ├── calibrate.R           # et_calibrate(), et_pit(), et_diagnose()
│   ├── plot.R                # et_plot_*() — eight ggplot2 visualizations
│   ├── data.R                # et_sim bundled dataset documentation
│   └── utils.R               # standardize(), unstandardize(), et_theme(), internal helpers
├── data/
│   └── et_sim.rda            # bundled simulated multi-cluster time-series dataset
├── tests/
│   └── testthat/
├── vignettes/
│   ├── introduction.Rmd      # Getting Started with ErrorTracer (full pipeline tutorial)
│   └── worked_example.Rmd    # A one-predictor worked example (companion to Box 1)
├── DESCRIPTION
└── NAMESPACE
```

---

## Dependencies

**Imports** (installed automatically):

| Package | Version | Purpose |
|---|---|---|
| `brms` | >= 2.20.0 | Bayesian regression via Stan |
| `ggplot2` | >= 3.4.0 | Visualization |
| `rlang` | >= 1.1.0 | Tidy evaluation (`.data` pronoun) |
| `stats`, `utils` | base R | Statistical functions and utilities |

**Suggests** (install separately as needed):

| Package | Version | Purpose |
|---|---|---|
| `glmnet` | >= 4.1.0 | Elastic net / lasso prior extraction |
| `ranger` | >= 0.15.0 | Random forest prior extraction (experimental) |
| `survival` | — | Kaplan–Meier estimation in `et_shelf_life_pool()` |
| `loo` | >= 2.6.0 | LOO cross-validation in `et_diagnose()` |
| `bayesplot` | >= 1.10.0 | Posterior predictive checks |
| `scales` | >= 1.3.0 | Percentage axis formatting in `et_plot_sensitivity()` |
| `tidyr`, `dplyr` | >= 1.3.0 / >= 1.1.0 | Data manipulation in vignettes |
| `testthat`, `covr` | >= 3.0.0 | Unit testing and coverage |
| `knitr`, `rmarkdown` | — | Vignette building |

**Backend** (required; install separately):
- `cmdstanr` (recommended) or `rstan`

---

## Caveats and Limitations

- **Small-n regime:** When `n < 10` observations, informed priors from `extract_priors()` dominate the posterior. This is the setting where the enet→prior pipeline provides the most value, but is also most sensitive to prior misspecification. Document prior choices carefully.
- **Intercept handling:** The first-stage intercept is excluded from the prior by default, and `intercept_prior_sd = NULL` lets `brms` use its data-aware default (`student_t(3, median(y), mad(y))`). Supplying a numeric `intercept_prior_sd` sets `Normal(0, intercept_prior_sd)`, which is only appropriate when the response is pre-centered.
- **Environmental uncertainty is inside the total:** `env_var` is estimated by perturbation and is a component *of* `total_var`, which since 1.4.0 is the measured `Var[posterior_predict]` (on the environmentally inclusive draws). Releases before 1.3.0 reported `env_var` outside the total, and 1.3.x defined the total as the component sum; decomposition tables from either should be recomputed rather than compared.
- **Predictions made before 1.4.1 are not reproducible:** `et_predict()` drew from the ambient random stream, so its intervals, PIT, calibration and shelf life depended on what had used that stream earlier in the session. It now seeds its draws from the model's own seed by default (`seed =` to override, `seed = NA` to opt out).
- **`ranger` priors are experimental and undirected:** Because random forests provide no signed coefficients or GLM-scale units, ranger-derived priors are centered at zero with heuristic importance-based SDs. Prefer `glmnet`/`lm`/`glm` for a principled prior.

---

## License

MIT License. See `LICENSE` for details.

---

## Use of AI assistance

Parts of this package were developed with the assistance of Anthropic's **Claude Opus**, used
as a coding and review tool. Its contributions include unit and regression tests, package and
function documentation, implementation work on several functions, and code audits that
identified and corrected defects in the variance decomposition and in the horizon-reporting
rules.

All AI-assisted output was specified, supervised, reviewed and accepted by the author, who
takes full responsibility for the correctness of the package. No AI system is an author or a
contributor to this software.

---

## Authors

- **Luis Javier Madrigal-Roca** (Creator, Maintainer) — madrigalrocalj@yahoo.com
- **John Kelly** (Author) — jkk@ku.edu
