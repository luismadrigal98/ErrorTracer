# R/utils.R — Internal utilities, logging, theming, numeric helpers

# ============================================================
# Internal logging (not exported)
# ============================================================

.et_log <- function(..., level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", ts, "] [ErrorTracer/", level, "] ", paste0(..., collapse = ""))
  message(msg)
}

.et_info  <- function(...) .et_log(..., level = "INFO")
.et_warn  <- function(...) .et_log(..., level = "WARN")
.et_error <- function(...) .et_log(..., level = "ERROR")

# ============================================================
# Numeric helpers
# ============================================================

#' Standardize a numeric vector (mean-centre, unit-variance)
#'
#' @param x Numeric vector.
#' @return Standardized numeric vector. Returns zeros if variance is zero.
#' @export
standardize <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s < 1e-12) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

#' Reverse standardization
#'
#' @param z Standardized values.
#' @param mu Original mean.
#' @param s  Original standard deviation.
#' @return Values on the original scale.
#' @export
unstandardize <- function(z, mu, s) {
  z * s + mu
}

# ============================================================
# ggplot2 theme
# ============================================================

#' Minimal ggplot2 theme for ErrorTracer plots
#'
#' @param base_size Base font size (default 12).
#' @return A \code{ggplot2} theme object.
#' @export
et_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 2),
      strip.text    = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom"
    )
}

# ============================================================
# Internal helpers for et_predict / decomposition
# ============================================================

# Normalise env_noise to a named list of per-predictor SDs.
# If scalar, replicate for every predictor in pred_names.
# If named list / named vector, validate names and fill missing with the
# default (0 means no perturbation for that predictor).
.resolve_env_noise <- function(env_noise, pred_names, newdata) {
  if (is.null(env_noise)) {
    noise_sds <- stats::setNames(rep(0, length(pred_names)), pred_names)
    return(noise_sds)
  }

  if (is.numeric(env_noise) && length(env_noise) == 1) {
    # Scalar: apply as a fraction of each predictor's SD in newdata
    noise_sds <- vapply(pred_names, function(p) {
      if (p %in% colnames(newdata)) {
        sd_p <- stats::sd(newdata[[p]], na.rm = TRUE)
        if (is.na(sd_p) || sd_p < 1e-12) sd_p <- 1
        env_noise * sd_p
      } else {
        0
      }
    }, numeric(1))
    return(noise_sds)
  }

  # Named list or named numeric vector
  env_noise <- unlist(env_noise)
  noise_sds <- stats::setNames(rep(0, length(pred_names)), pred_names)
  shared <- intersect(names(env_noise), pred_names)
  noise_sds[shared] <- env_noise[shared]
  if (length(shared) < length(names(env_noise))) {
    unknown <- setdiff(names(env_noise), pred_names)
    .et_warn("env_noise contains predictor(s) not in model: ",
             paste(unknown, collapse = ", "), " — ignored")
  }
  noise_sds
}

# Extract the names of fixed-effect predictors (excludes Intercept) from a
# brmsfit object.
.brms_pred_names <- function(fit) {
  fe <- rownames(brms::fixef(fit))
  setdiff(fe, "Intercept")
}

# Extract posterior draws matrix from a brmsfit (rows = draws, cols = params).
# Optionally thin to at most max_draws rows.
.brms_draws_matrix <- function(fit, max_draws = NULL) {
  mat <- as.matrix(fit)
  if (!is.null(max_draws) && nrow(mat) > max_draws) {
    mat <- mat[seq_len(max_draws), , drop = FALSE]
  }
  mat
}
