# R/decompose.R — decompose_uncertainty(): three-way variance decomposition

#' Extract or recompute uncertainty decomposition
#'
#' Returns a \code{data.frame} with the three-way uncertainty decomposition
#' stored inside an \code{et_prediction} object:
#' \describe{
#'   \item{param_var}{Variance of the posterior linear predictor — captures
#'     uncertainty in fitted regression coefficients.}
#'   \item{env_var}{Additional variance arising from measurement or
#'     prediction uncertainty in the predictor values (estimated via
#'     perturbation in \code{\link{et_predict}}).  Zero when
#'     \code{env_noise = NULL}.}
#'   \item{residual_var}{Posterior mean of \eqn{\sigma^2} — biological
#'     process noise, unmeasured drivers, and drift. Using the mean (not
#'     median) ensures \code{param_var + residual_var} \eqn{\approx}{~} \code{total_var}
#'     by the law of total variance.}
#'   \item{total_var}{Variance of the full posterior predictive draws
#'     (parameter + residual; note that \code{env_var} is an additive
#'     component measured separately from the perturbation step).}
#' }
#'
#' All variance components are guaranteed non-negative.
#'
#' @param predictions An \code{et_prediction} object from
#'   \code{\link{et_predict}}, or an \code{et_prediction_list} (grouped).
#' @param ... Unused.
#' @return A \code{data.frame} with columns
#'   \code{obs_id, param_var, env_var, residual_var, total_var}.
#'   For grouped predictions, a \code{group} column is prepended.
#' @examples
#' \donttest{
#' set.seed(1)
#' df  <- data.frame(y = rnorm(20), x1 = rnorm(20))
#' fit <- et_fit(y ~ x1, data = df,
#'               chains = 1, iter = 500, warmup = 250,
#'               cores = 1, refresh = 0)
#' new_df <- data.frame(x1 = rnorm(5))
#' pred   <- et_predict(fit, newdata = new_df,
#'                      env_noise = list(x1 = 0.2),
#'                      n_draws = 200, n_perturb = 50)
#' decomp <- decompose_uncertainty(pred)
#' head(decomp)
#' }
#' @export
decompose_uncertainty <- function(predictions, ...) {
  UseMethod("decompose_uncertainty")
}

#' @export
decompose_uncertainty.et_prediction <- function(predictions, ...) {
  predictions$decomposition
}

#' @export
decompose_uncertainty.et_prediction_list <- function(predictions, ...) {
  parts <- lapply(names(predictions$predictions), function(g) {
    pred <- predictions$predictions[[g]]
    if (is.null(pred)) return(NULL)
    d <- pred$decomposition
    cbind(data.frame(group = g, stringsAsFactors = FALSE), d)
  })
  do.call(rbind, Filter(Negate(is.null), parts))
}

#' @export
decompose_uncertainty.default <- function(predictions, ...) {
  stop("decompose_uncertainty() expects an et_prediction or et_prediction_list object.")
}
