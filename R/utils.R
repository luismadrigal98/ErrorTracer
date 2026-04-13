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

# Normalise env_noise to a named list of per-observation SD vectors.
# Each element of the returned list has length n_obs, so every row of
# newdata can have its own noise SD (enabling time-varying uncertainty).
#
# Accepted forms for env_noise:
#   NULL                  — zero noise for all predictors
#   scalar numeric        — fraction of each predictor's SD in newdata
#   named scalar list/vec — fixed SD per predictor (constant across obs)
#   named vector list     — per-row SDs; each entry must be length 1 or n_obs
.resolve_env_noise <- function(env_noise, pred_names, newdata) {
  n_obs <- nrow(newdata)
  zeros <- stats::setNames(lapply(pred_names, function(p) rep(0, n_obs)),
                           pred_names)

  if (is.null(env_noise)) return(zeros)

  if (is.numeric(env_noise) && length(env_noise) == 1) {
    # Scalar fraction: scale by each predictor's empirical SD
    result <- lapply(pred_names, function(p) {
      sd_p <- if (p %in% colnames(newdata)) {
        s <- stats::sd(newdata[[p]], na.rm = TRUE)
        if (is.na(s) || s < 1e-12) 1 else s
      } else 1
      rep(env_noise * sd_p, n_obs)
    })
    return(stats::setNames(result, pred_names))
  }

  # Named list or named numeric vector (possibly with per-row vectors)
  env_noise_list <- if (is.list(env_noise)) env_noise else as.list(env_noise)

  unknown <- setdiff(names(env_noise_list), pred_names)
  if (length(unknown) > 0) {
    .et_warn("env_noise contains predictor(s) not in model: ",
             paste(unknown, collapse = ", "), " — ignored")
  }

  result <- zeros
  shared <- intersect(names(env_noise_list), pred_names)

  for (p in shared) {
    v <- as.numeric(env_noise_list[[p]])
    if (length(v) == 1L) {
      result[[p]] <- rep(v, n_obs)
    } else if (length(v) == n_obs) {
      result[[p]] <- v
    } else {
      stop(
        "env_noise[[\"", p, "\"]] has length ", length(v),
        " but newdata has ", n_obs, " row(s). ",
        "Supply a scalar (constant noise) or a vector of length n_obs ",
        "(time-varying noise)."
      )
    }
  }

  result
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
