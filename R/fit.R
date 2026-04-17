# R/fit.R — et_fit(): Bayesian model fitting via brms

# ============================================================
# Main fitting function
# ============================================================

#' Fit a Bayesian regression model with informed priors
#'
#' Wraps \code{brms::brm()} and attaches the prior specification,
#' training data reference, and configuration for downstream uncertainty
#' decomposition.  Pass \code{priors} from \code{\link{extract_priors}} to
#' use regularized-model coefficients as prior means; omit it for
#' default (weakly informative) priors.
#'
#' @param formula An R formula, e.g.\ \code{response ~ .} or
#'   \code{y ~ x1 + x2}.
#' @param data A \code{data.frame} with all predictors and the response.
#' @param priors An \code{et_prior_spec} object from
#'   \code{\link{extract_priors}}, or a \code{brmsprior} object, or
#'   \code{NULL} for brms defaults.
#' @param chains Integer.  Number of MCMC chains (default 4).
#' @param iter Integer.  Total iterations per chain, including warmup
#'   (default 2000).
#' @param warmup Integer.  Warmup iterations per chain (default
#'   \code{floor(iter / 2)}).
#' @param cores Integer.  Parallel cores (default
#'   \code{min(chains, parallel::detectCores())}).
#' @param seed Integer.  Random seed for reproducibility (default 42).
#' @param adapt_delta Numeric.  Target acceptance probability for HMC
#'   (default 0.95).
#' @param max_treedepth Integer.  Maximum tree depth (default 12).
#' @param grouping Character.  Name of a column in \code{data} to use for
#'   grouping.  If non-\code{NULL}, one model is fitted per unique group
#'   value and an \code{et_model_list} is returned.
#' @param silent Integer passed to \code{brms::brm()} (default 2, no Stan
#'   output).
#' @param ... Additional arguments passed to \code{brms::brm()}.
#' @return An \code{et_model} object (or an \code{et_model_list} if
#'   \code{grouping} is specified).
#' @examples
#' \dontrun{
#' ps  <- extract_priors(lm(mpg ~ wt + hp, data = mtcars))
#' fit <- et_fit(mpg ~ wt + hp, data = mtcars, priors = ps, chains = 2, iter = 1000)
#' }
#' @export
et_fit <- function(formula,
                   data,
                   priors = NULL,
                   chains = 4L,
                   iter = 2000L,
                   warmup = floor(iter / 2),
                   cores = min(chains, parallel::detectCores()),
                   seed = 42L,
                   adapt_delta = 0.95,
                   max_treedepth = 12L,
                   grouping = NULL,
                   silent = 2L,
                   ...) {

  if (!is.null(grouping)) {
    return(.et_fit_grouped(
      formula = formula, data = data, priors = priors,
      chains = chains, iter = iter, warmup = warmup, cores = cores,
      seed = seed, adapt_delta = adapt_delta, max_treedepth = max_treedepth,
      grouping = grouping, silent = silent, ...
    ))
  }

  .et_fit_single(
    formula = formula, data = data, priors = priors,
    chains = chains, iter = iter, warmup = warmup, cores = cores,
    seed = seed, adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    silent = silent, ...
  )
}

# ============================================================
# Internal: fit a single model
# ============================================================

.et_fit_single <- function(formula, data, priors, chains, iter, warmup,
                            cores, seed, adapt_delta, max_treedepth,
                            silent, ...) {

  brms_prior <- if (inherits(priors, "et_prior_spec")) priors$prior
                else priors  # brmsprior or NULL

  n_pred <- length(attr(stats::terms(formula, data = data), "term.labels"))
  .et_info("Fitting Bayesian model: ", deparse(formula),
           " (", nrow(data), " obs, ~", n_pred, " predictors)")

  fit <- brms::brm(
    formula = brms::bf(formula),
    data = data,
    prior = brms_prior,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    silent  = silent,
    refresh = 0,
    ...
  )

  structure(
    list(
      fit = fit,
      prior_spec = if (inherits(priors, "et_prior_spec")) priors else NULL,
      formula = formula,
      data = data,
      config = list(
        chains = chains, iter = iter, warmup = warmup,
        cores = cores,  seed = seed,
        adapt_delta = adapt_delta, max_treedepth = max_treedepth
      )
    ),
    class = "et_model"
  )
}

# ============================================================
# Internal: fit one model per group
# ============================================================

.et_fit_grouped <- function(formula, data, priors, chains, iter, warmup,
                             cores, seed, adapt_delta, max_treedepth,
                             grouping, silent, ...) {

  if (!grouping %in% colnames(data)) {
    stop("grouping column '", grouping, "' not found in data.")
  }

  groups <- unique(data[[grouping]])
  models <- vector("list", length(groups))
  names(models) <- as.character(groups)

  for (g in groups) {
    gname <- as.character(g)
    .et_info("Fitting model for group: ", gname)
    sub_data <- data[data[[grouping]] == g, , drop = FALSE]

    # Use group-specific prior if priors is a named list
    group_prior <- if (is.list(priors) && !inherits(priors, "et_prior_spec")) {
      priors[[gname]]
    } else {
      priors
    }

    models[[gname]] <- tryCatch(
      .et_fit_single(
        formula = formula, data = sub_data, priors = group_prior,
        chains = chains, iter = iter, warmup = warmup, cores = cores,
        seed = seed, adapt_delta = adapt_delta, max_treedepth = max_treedepth,
        silent = silent, ...
      ),
      error = function(e) {
        .et_error("Failed to fit model for group ", gname, ": ", e$message)
        NULL
      }
    )
  }

  structure(
    list(
      models   = models,
      grouping = grouping,
      formula  = formula
    ),
    class = "et_model_list"
  )
}

# ============================================================
# S3 methods for et_model
# ============================================================

#' @export
print.et_model <- function(x, ...) {
  cat("ErrorTracer model (et_model)\n")
  cat("  Formula :", deparse(x$formula), "\n")
  cat("  n obs   :", nrow(x$data), "\n")
  cat("  Chains  :", x$config$chains,
      "  Iter:", x$config$iter,
      "  Warmup:", x$config$warmup, "\n")
  if (!is.null(x$prior_spec)) {
    cat("  Priors  : informed (", x$prior_spec$method, ", ",
        length(x$prior_spec$pred_names), " predictors)\n", sep = "")
  } else {
    cat("  Priors  : brms defaults\n")
  }
  rhat_max <- tryCatch(max(brms::rhat(x$fit), na.rm = TRUE), error = function(e) NA)
  cat("  Rhat max:", if (is.na(rhat_max)) "NA" else round(rhat_max, 3), "\n")
  invisible(x)
}

#' @export
summary.et_model <- function(object, ...) {
  cat("=== ErrorTracer model summary ===\n\n")
  print(object)
  cat("\n--- Fixed effects ---\n")
  print(brms::fixef(object$fit))
  invisible(object)
}

# ============================================================
# S3 methods for et_model_list
# ============================================================

#' @export
print.et_model_list <- function(x, ...) {
  cat("ErrorTracer grouped model list (et_model_list)\n")
  cat("  Grouping :", x$grouping, "\n")
  cat("  Formula  :", deparse(x$formula), "\n")
  cat("  Groups   :", length(x$models), "\n")
  fitted <- sum(!vapply(x$models, is.null, logical(1)))
  cat("  Fitted   :", fitted, "/", length(x$models), "\n")
  invisible(x)
}

#' @export
summary.et_model_list <- function(object, ...) {
  print(object)
  cat("\n--- Per-group Rhat max ---\n")
  for (nm in names(object$models)) {
    m <- object$models[[nm]]
    if (is.null(m)) {
      cat("  ", nm, ": FAILED\n")
    } else {
      rhat_max <- tryCatch(max(brms::rhat(m$fit), na.rm = TRUE), error = function(e) NA)
      cat(sprintf("  %-20s  Rhat max = %.3f\n", nm,
                  if (is.na(rhat_max)) NA else rhat_max))
    }
  }
  invisible(object)
}
