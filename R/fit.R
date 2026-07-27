# R/fit.R — et_fit(): Bayesian model fitting via brms

# ******************************************************************************
# Main fitting function
# ______________________________________________________________________________

#' Fit a Bayesian regression model with informed priors
#'
#' Wraps \code{brms::brm()} and attaches the prior specification,
#' training data reference, and configuration for downstream uncertainty
#' decomposition.  Pass \code{priors} from \code{\link{extract_priors}} to
#' use regularized-model coefficients to set the prior scale (and, with
#' \code{shrinkage = "estimate"}, the prior mean); omit it for default
#' (weakly informative) priors.
#'
#' @section Scope and limitations:
#' \code{et_fit()} is deliberately a \strong{thin wrapper around \code{brms}}:
#' the formula, family, and any autocorrelation term (\code{ar()}/\code{ma()}/
#' \code{arma()}/...) are passed straight through, so ErrorTracer inherits
#' \code{brms}'s modelling scope and its Stan back-end. The added value is the
#' prior pipeline, the structured uncertainty decomposition, forecast skill /
#' shelf life, and calibration diagnostics --- not new model classes.
#'
#' Two boundaries to keep in mind:
#' \itemize{
#'   \item \strong{Hierarchical / random-effects models} fit through the
#'     formula (\code{y ~ x + (x | group)}), but the current decomposition does
#'     \emph{not} yet split out a group-level variance component, and it does
#'     not distinguish in-sample prediction (using group-specific parameters)
#'     from out-of-sample prediction (integrating over the group level). Treat
#'     the decomposition of hierarchical fits as provisional until a dedicated
#'     group-variance term is added. The \code{grouping} argument here is a
#'     \emph{separate} device: it fits one independent model per group (no
#'     pooling), not a single multilevel model.
#'   \item \strong{Convergence is your responsibility.} \code{et_fit()} emits a
#'     warning on high Rhat or divergent transitions, but you should inspect
#'     \code{\link{et_diagnose}} before predicting; the defaults
#'     (\code{iter = 2000}, 4 chains) are a starting point, not a guarantee.
#' }
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
#' @param eiv Optional errors-in-variables specification.  A named list /
#'   vector mapping predictor names to either a scalar SD or a vector of
#'   per-row SDs (length \code{nrow(data)}).  For each entry, the formula
#'   term for that predictor is rewritten as \code{brms::me(pred, se_pred)}
#'   (an auxiliary \code{se_<pred>} column is appended to \code{data}), so
#'   the posterior reflects measurement error in the predictor as well as
#'   coefficient uncertainty.  The beta posteriors widen accordingly, which
#'   partially absorbs what ErrorTracer's downstream \code{env_var}
#'   component would otherwise report.  When \code{eiv} is supplied together
#'   with an \code{et_prior_spec} from \code{\link{extract_priors}}, the
#'   informed priors are \emph{dropped} because they target \code{class = "b"}
#'   terms and \code{me()} terms live under \code{class = "bsp"}; brms
#'   defaults are used instead (and a warning is logged).
#' @param ar_prior Character. Prior on autocorrelation parameters
#'   (\code{class = "ar"} / \code{"ma"}) when the formula carries an
#'   \code{ar()} / \code{ma()} / \code{arma()} term. \code{brms}'s own default
#'   is \emph{flat and unbounded}, which on a short series lets the posterior
#'   straddle the unit root; the k-step forecast variance
#'   \eqn{\sigma^2 (1 - \phi^{2k}) / (1 - \phi^2)} then diverges, and
#'   \code{temporal_var} and \code{\link{shelf_life}} become uninterpretable.
#'   Options:
#'   \describe{
#'     \item{\code{"weakly_informative"}}{(default) \code{normal(0, 0.5)} — about
#'       95\% of the prior mass lies inside \eqn{(-1, 1)}, but the parameter is
#'       not hard-bounded, so a genuinely non-stationary series can still say so.}
#'     \item{\code{"stationary"}}{The same prior truncated to \eqn{(-1, 1)}. For
#'       AR(1) this is exactly the stationarity region; for higher orders it is
#'       necessary but not sufficient, and a warning says so.}
#'     \item{\code{"flat"}}{\code{brms}'s unbounded default.}
#'   }
#'   A prior the caller supplies for these classes is never overridden.
#'   Regardless of the setting, \code{et_fit()} warns after fitting when more
#'   than 5\% of the posterior sits at \eqn{|\phi| \ge 1}.
#' @param silent Integer passed to \code{brms::brm()} (default 2, no Stan
#'   output).
#' @param ... Additional arguments passed to \code{brms::brm()}.
#' @return An \code{et_model} object (or an \code{et_model_list} if
#'   \code{grouping} is specified).
#' @examples
#' \donttest{
#' set.seed(1)
#' df  <- data.frame(y = rnorm(20), x1 = rnorm(20), x2 = rnorm(20))
#' ps  <- extract_priors(lm(y ~ x1 + x2, data = df))
#' fit <- et_fit(y ~ x1 + x2, data = df, priors = ps,
#'               chains = 1, iter = 500, warmup = 250,
#'               cores = 1, refresh = 0)
#' print(fit)
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
                   eiv = NULL,
                   ar_prior = c("weakly_informative", "stationary", "flat"),
                   silent = 2L,
                   ...) {

  ar_prior <- match.arg(ar_prior)

  # Rewrite formula + augment data for measurement-error predictors if eiv
  # was supplied. Strips any et_prior_spec priors because their (class="b",
  # coef=pred_name) entries don't apply to me() (class="bsp") terms.
  eiv_spec <- NULL
  if (!is.null(eiv)) {
    rewrite  <- .apply_eiv(formula, data, eiv)
    formula  <- rewrite$formula
    data     <- rewrite$data
    eiv_spec <- rewrite$eiv_spec
    # brms parses me() via eval() against the search path, so brms must be
    # attached (not just imported) for me() to resolve. Attach it once.
    if (!"package:brms" %in% search()) {
      suppressPackageStartupMessages(attachNamespace("brms"))
    }
    if (inherits(priors, "et_prior_spec")) {
      .et_warn("eiv specified; dropping informed priors (they target ",
               "class='b' terms but me() coefficients are class='bsp'). ",
               "Using brms defaults.")
      priors <- NULL
    }
  }

  if (!is.null(grouping)) {
    result <- .et_fit_grouped(
      formula = formula, data = data, priors = priors,
      chains = chains, iter = iter, warmup = warmup, cores = cores,
      seed = seed, adapt_delta = adapt_delta, max_treedepth = max_treedepth,
      grouping = grouping, ar_prior = ar_prior, silent = silent, ...
    )
    if (!is.null(eiv_spec)) result$eiv_spec <- eiv_spec
    return(result)
  }

  result <- .et_fit_single(
    formula = formula, data = data, priors = priors,
    chains = chains, iter = iter, warmup = warmup, cores = cores,
    seed = seed, adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    ar_prior = ar_prior, silent = silent, ...
  )
  if (!is.null(eiv_spec)) result$eiv_spec <- eiv_spec
  result
}

# ******************************************************************************
# Internal: rewrite formula and augment data for errors-in-variables
# ______________________________________________________________________________

# Given (formula, data, eiv), return a list(formula = new, data = new).
# For each named predictor in eiv:
#   * append a column se_<pred> to data (scalar recycled or vector of length
#     nrow(data)),
#   * substitute pred -> me(pred, se_<pred>) in the formula RHS.
# The substitution preserves everything else in the formula (intercept,
# interactions, factor contrasts). Matches only on whole-token predictor
# names to avoid partial-string collisions.
.apply_eiv <- function(formula, data, eiv) {
  if (!is.list(eiv) && !is.numeric(eiv)) {
    stop("eiv must be a named list or named numeric vector.")
  }
  if (is.null(names(eiv)) || any(!nzchar(names(eiv)))) {
    stop("eiv must have names matching predictor columns.")
  }
  eiv <- as.list(eiv)

  term_labels <- attr(stats::terms(formula, data = data), "term.labels")
  missing_preds <- setdiff(names(eiv), colnames(data))
  if (length(missing_preds)) {
    stop("eiv references column(s) not in data: ",
         paste(missing_preds, collapse = ", "))
  }

  n_obs <- nrow(data)
  for (p in names(eiv)) {
    v <- as.numeric(eiv[[p]])
    if (length(v) == 1L) v <- rep(v, n_obs)
    if (length(v) != n_obs) {
      stop("eiv[['", p, "']] has length ", length(v),
           " but data has ", n_obs, " row(s).")
    }
    if (any(is.na(v) | v < 0)) {
      stop("eiv[['", p, "']] must be non-negative and finite.")
    }
    data[[paste0("se_", p)]] <- v
  }

  rhs <- deparse(formula[[3L]], width.cutoff = 500L)
  rhs <- paste(rhs, collapse = " ")
  for (p in names(eiv)) {
    pattern <- paste0("(?<![A-Za-z0-9_.])", p, "(?![A-Za-z0-9_.])")
    replacement <- paste0("me(", p, ", se_", p, ")")
    rhs <- gsub(pattern, replacement, rhs, perl = TRUE)
  }
  lhs <- deparse(formula[[2L]], width.cutoff = 500L)
  # brms evaluates `me()` during formula parsing via eval() on the formula's
  # own environment, so we bind `me` (and `mi`, for completeness) directly
  # in that environment rather than relying on brms being attached.
  env <- new.env(parent = environment(formula) %||% parent.frame())
  env$me <- brms::me
  if (exists("mi", envir = asNamespace("brms"))) env$mi <- brms::mi
  new_formula <- stats::as.formula(paste(lhs, "~", rhs), env = env)

  list(formula = new_formula, data = data, eiv_spec = eiv)
}

# ******************************************************************************
# Internal: fit a single model
# ______________________________________________________________________________

.et_fit_single <- function(formula, data, priors, chains, iter, warmup,
                            cores, seed, adapt_delta, max_treedepth,
                            silent, ar_prior = "weakly_informative", ...) {

  brms_prior <- if (inherits(priors, "et_prior_spec")) priors$prior
                else priors  # brmsprior or NULL

  # Autocorrelation parameters: brms's default prior on class "ar" / "ma" is
  # FLAT AND UNBOUNDED, so nothing keeps the process inside the stationary
  # region. On a short series that routinely yields posteriors straddling the
  # unit root, whose k-step forecast variance then explodes. Attach a
  # weakly-informative default instead (see .ar_prior_spec).
  ac_prior <- .ar_prior_spec(formula, data, ar_prior, existing = brms_prior)
  if (!is.null(ac_prior)) {
    brms_prior <- if (is.null(brms_prior)) ac_prior else c(brms_prior, ac_prior)
  }

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

  # Convergence nudge: a quick Rhat / divergence check so the user is warned
  # BEFORE predicting (Stage 3). Full diagnostics live in et_diagnose().
  .et_check_convergence(fit)
  # Stationarity nudge: a near-unit-root autocorrelation posterior makes the
  # far-lead predictive variance heavy-tailed, which downstream shows up as a
  # dominant temporal_var and an exploding shelf-life ratio. Surface it here,
  # where the user can still act on it, rather than at decomposition time.
  .et_check_stationarity(fit)

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

# ******************************************************************************
# Internal: autocorrelation-parameter priors and stationarity diagnostic
# ______________________________________________________________________________

# brms puts a FLAT, UNBOUNDED prior on class "ar" / "ma"
# (get_prior() reports "(flat)" with no lb/ub). For a long, well-identified
# series that is harmless; on the short series typical of ecological panels it
# lets the posterior straddle |phi| >= 1, and the k-step-ahead forecast variance
# sigma^2 (1 - phi^(2k)) / (1 - phi^2) then diverges. The downstream symptoms are
# a heavy-tailed predictive distribution, a temporal_var that swamps every other
# channel, and a shelf-life ratio that explodes.
#
# Defaults:
#   "weakly_informative" (default) -- normal(0, 0.5), which places ~95% of the
#       prior mass inside (-1, 1) but does NOT hard-bound the parameter, so a
#       genuinely non-stationary series can still say so.
#   "stationary" -- the same prior truncated to (-1, 1). For AR(1) this IS the
#       stationarity region. For p > 1 it is necessary but NOT sufficient, so we
#       say so rather than implying a guarantee.
#   "flat" -- brms's default; preserves the pre-existing behaviour.
#
# Never overrides a class the caller already specified.
.ar_prior_spec <- function(formula, data, ar_prior = "weakly_informative",
                           existing = NULL) {
  if (identical(ar_prior, "flat")) return(NULL)
  if (!.formula_has_autocor(formula)) return(NULL)

  classes <- tryCatch({
    gp <- brms::get_prior(brms::bf(formula), data = data)
    intersect(unique(as.character(gp$class)), c("ar", "ma"))
  }, error = function(e) character(0))
  if (!length(classes)) return(NULL)

  if (!is.null(existing) && !is.null(existing$class)) {
    classes <- setdiff(classes, as.character(existing$class))
  }
  if (!length(classes)) return(NULL)

  stationary <- identical(ar_prior, "stationary")

  if (stationary) {
    # Only claim a stationarity guarantee where one actually holds.
    p_order <- tryCatch({
      ac  <- brms::brmsterms(formula)$dpars$mu$ac
      ord <- suppressWarnings(as.numeric(c(ac$p, ac$q)))
      ord <- ord[is.finite(ord)]
      if (length(ord)) max(ord) else NA_real_
    }, error = function(e) NA_real_)
    if (is.finite(p_order) && p_order > 1) {
      .et_warn("ar_prior = 'stationary' bounds each autocorrelation ",
               "coefficient to (-1, 1). For order ", p_order,
               " that is necessary but NOT sufficient for stationarity; the ",
               "stationary region is not the hypercube. Check the fitted ",
               "coefficients before trusting long-lead forecasts.")
    }
  }

  parts <- lapply(classes, function(cl) {
    if (stationary) {
      brms::set_prior("normal(0, 0.5)", class = cl, lb = -1, ub = 1)
    } else {
      brms::set_prior("normal(0, 0.5)", class = cl)
    }
  })
  do.call(c, parts)
}

# Post-fit stationarity nudge. For AR(1), P(|phi| >= 1) is exactly the posterior
# probability of a non-stationary process; for higher orders it is a coarse
# proxy, which the message reflects. Warns, never stops -- a random-walk residual
# can be the correct answer, but the user must know it is what they have.
.et_check_stationarity <- function(fit) {
  dm <- tryCatch(as.matrix(fit), error = function(e) NULL)
  if (is.null(dm)) return(invisible(NULL))
  ac_cols <- grep("^(ar|ma)\\[", colnames(dm), value = TRUE)
  if (!length(ac_cols)) return(invisible(NULL))

  for (cl in ac_cols) {
    p_ge1 <- mean(abs(dm[, cl]) >= 1, na.rm = TRUE)
    if (is.finite(p_ge1) && p_ge1 > 0.05) {
      .et_warn("Autocorrelation parameter ", cl, " has ",
               round(100 * p_ge1), "% of its posterior at |value| >= 1 ",
               "(posterior mean ", round(mean(dm[, cl], na.rm = TRUE), 3),
               "). The residual process is at or beyond the unit root, so the ",
               "k-step forecast variance grows without bound and long-lead ",
               "intervals, temporal_var and shelf_life() should be read as ",
               "'no usable horizon' rather than as calibrated numbers. ",
               "Consider ar_prior = 'stationary', a longer training series, ",
               "or differencing the response.")
    }
  }
  invisible(NULL)
}

# Lightweight post-fit convergence nudge. Warns (never stops) on high Rhat or
# divergent transitions and points the user to et_diagnose() before Stage 3.
.et_check_convergence <- function(fit) {
  info <- tryCatch({
    rh   <- suppressWarnings(brms::rhat(fit))
    nuts <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
    divs <- if (!is.null(nuts))
      sum(nuts$Value[nuts$Parameter == "divergent__"], na.rm = TRUE) else 0
    list(max_rh = max(rh, na.rm = TRUE), divs = divs)
  }, error = function(e) NULL)
  if (is.null(info)) return(invisible())

  msgs <- character(0)
  if (is.finite(info$max_rh) && info$max_rh > 1.05) {
    msgs <- c(msgs, sprintf("max Rhat = %.3f (> 1.05)", info$max_rh))
  }
  if (info$divs > 0) {
    msgs <- c(msgs, sprintf("%d divergent transition(s)", info$divs))
  }
  if (length(msgs)) {
    .et_warn("Convergence warning: ", paste(msgs, collapse = "; "),
             ". Inspect with et_diagnose() and consider more iterations or a ",
             "higher adapt_delta BEFORE predicting (Stage 3).")
  }
  invisible()
}

# ******************************************************************************
# Internal: fit one model per group
# ______________________________________________________________________________

.et_fit_grouped <- function(formula, data, priors, chains, iter, warmup,
                             cores, seed, adapt_delta, max_treedepth,
                             grouping, silent,
                             ar_prior = "weakly_informative", ...) {

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
        ar_prior = ar_prior, silent = silent, ...
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

# ******************************************************************************
# S3 methods for et_model
# ______________________________________________________________________________

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

# ******************************************************************************
# S3 methods for et_model_list
# ______________________________________________________________________________

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
