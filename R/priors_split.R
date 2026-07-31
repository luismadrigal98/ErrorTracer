# R/priors_split.R -- data-splitting prior extraction (removes double use)

#' Extract priors from data held out from the likelihood
#'
#' The standard workflow fits a regularized model and the Bayesian model on the
#' \emph{same} rows.  Even with \code{shrinkage = "zero"}, which never reuses a
#' point estimate as a prior mean, two things still cross over: the prior
#' \emph{scale} is set from those rows, and --- more consequentially --- when the
#' first-stage model performs variable \emph{selection}, the surviving predictor
#' set is chosen using the same data that then forms the likelihood.  Reported
#' coverage is therefore conditional on that selection being correct, which is
#' easy to satisfy in a simulation with orthogonal nuisance predictors and much
#' less so in the correlated designs typical of real covariates.
#'
#' \code{et_priors_split()} removes both by construction: the first-stage model
#' is fitted on a random subset, and the rows returned for the Bayesian fit are
#' its complement.  Selection and likelihood then share no observations, so the
#' posterior needs no post-selection caveat.
#'
#' @section The trade-off, stated plainly:
#' Splitting costs sample size.  With \code{prop = 0.5} the likelihood sees half
#' the data, so the posterior is genuinely wider --- that width is honest rather
#' than pessimistic, but at small \eqn{n} the loss can dominate whatever the
#' prior contributes.  Splitting is worth it when selection stability is in
#' doubt (many correlated candidate predictors) and not worth it when the
#' predictor set is fixed a priori, where there is no selection to protect
#' against and the full sample should go to the likelihood.
#'
#' @param data A \code{data.frame} containing both the first-stage predictors
#'   and the response.
#' @param prior_fun A function taking a \code{data.frame} and returning a fitted
#'   model that \code{\link{extract_priors}} accepts (e.g.\ an
#'   \code{lm}/\code{glm}/\code{cv.glmnet} fit).
#' @param prop Proportion of rows used for the \strong{prior} stage (default
#'   \code{0.5}).  The remainder is returned for the likelihood.
#' @param seed Optional integer seed for the split, for reproducibility.
#' @param strata Optional column name in \code{data}.  When supplied the split
#'   is stratified within its levels, so each level is represented in both
#'   halves --- important for grouped or seasonal designs where a naive split
#'   can leave a level entirely absent from one side.
#' @param ... Passed to \code{\link{extract_priors}} (e.g.\ \code{multiplier},
#'   \code{min_sd}, \code{shrinkage}).
#'
#' @return A list with:
#'   \describe{
#'     \item{priors}{The \code{et_priors} object from the prior-stage rows.}
#'     \item{fit_data}{Rows to pass to \code{\link{et_fit}} --- disjoint from
#'       the rows used to build the priors.}
#'     \item{prior_data}{Rows used for the first stage, returned so the split is
#'       auditable rather than hidden.}
#'     \item{prop, n_prior, n_fit}{The realised split.}
#'   }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
#' d$y <- 1 + 2 * d$x1 + rnorm(80)
#' sp  <- et_priors_split(d, function(df) lm(y ~ x1 + x2, data = df), seed = 1)
#' fit <- et_fit(y ~ x1 + x2, data = sp$fit_data, priors = sp$priors,
#'               chains = 1, iter = 500, refresh = 0)
#' }
#' @seealso \code{\link{extract_priors}} for the single-sample workflow.
#' @export
et_priors_split <- function(data, prior_fun, prop = 0.5, seed = NULL,
                            strata = NULL, ...) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.function(prior_fun)) stop("`prior_fun` must be a function.")
  if (!is.numeric(prop) || prop <= 0 || prop >= 1) {
    stop("`prop` must be strictly between 0 and 1.")
  }
  n <- nrow(data)
  if (n < 4L) stop("Need at least 4 rows to split.")

  if (!is.null(seed)) {
    old <- .Random.seed_state()
    on.exit(.Random.seed_restore(old), add = TRUE)
    set.seed(seed)
  }

  idx_prior <- if (is.null(strata)) {
    sort(sample.int(n, size = max(2L, floor(prop * n))))
  } else {
    if (!strata %in% colnames(data)) {
      stop("`strata` column '", strata, "' not found in data.")
    }
    lev <- split(seq_len(n), data[[strata]])
    sort(unlist(lapply(lev, function(ix) {
      if (length(ix) < 2L) return(ix)   # too small to split; keep on prior side
      sample(ix, size = max(1L, floor(prop * length(ix))))
    }), use.names = FALSE))
  }

  prior_data <- data[idx_prior, , drop = FALSE]
  fit_data   <- data[-idx_prior, , drop = FALSE]

  if (nrow(fit_data) < 2L) {
    stop("Split leaves ", nrow(fit_data), " row(s) for the likelihood; ",
         "lower `prop`.")
  }

  stage1 <- prior_fun(prior_data)
  priors <- extract_priors(stage1, ...)

  .et_info("Prior/likelihood split: ", nrow(prior_data), " row(s) for the ",
           "first stage, ", nrow(fit_data), " disjoint row(s) for the ",
           "likelihood. Selection and likelihood share no observations, so no ",
           "post-selection caveat applies to the resulting posterior.")

  list(
    priors     = priors,
    fit_data   = fit_data,
    prior_data = prior_data,
    prop       = nrow(prior_data) / n,
    n_prior    = nrow(prior_data),
    n_fit      = nrow(fit_data)
  )
}

# Save/restore the RNG state so a user-supplied seed does not leak into the
# caller's stream (set.seed() inside a function otherwise silently reseeds the
# session, which makes downstream "reproducible" results depend on whether this
# function was called).
.Random.seed_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

.Random.seed_restore <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
}
