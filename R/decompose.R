# R/decompose.R — decompose_uncertainty(): variance decomposition
# Reports a three-way split (parameter / environmental / residual) for iid
# models, and a four-way split adding a temporal-autocorrelation component
# when the model formula contains an ar() / ma() / arma() / cosy() / unstr()
# / sar() / car() term.
#
# total_var is the MEASURED variance of the env-inclusive posterior predictive
# draws, not the sum of the channels. See .decompose_from_arrays() in predict.R.

#' Extract or recompute uncertainty decomposition
#'
#' Returns a \code{data.frame} with the uncertainty decomposition stored
#' inside an \code{et_prediction} object:
#' \describe{
#'   \item{param_var}{Variance of the posterior linear predictor — captures
#'     uncertainty in fitted regression coefficients.}
#'   \item{env_var}{Variance arising from measurement or prediction
#'     uncertainty in the predictor values (estimated via perturbation in
#'     \code{\link{et_predict}}).  A genuine sub-share of \code{total_var}
#'     (see below), not an add-on. Zero when \code{env_noise = NULL}.}
#'   \item{residual_var}{Posterior mean of \eqn{\sigma^2} (or its
#'     family-specific analogue) — biological process noise, unmeasured
#'     drivers, and drift. For autocorrelation models this is the
#'     \emph{innovation} variance, not the stationary marginal variance;
#'     the autocorrelated accumulation is reported separately in
#'     \code{temporal_var}.  For a non-identity link it is evaluated at the
#'     \emph{perturbed} mean, \eqn{E_s[V(\mu^{\mathrm{env}})]}, because that is
#'     the mean the predictive draws are generated at; it therefore depends on
#'     \code{env_noise} through the variance function.  For Gaussian,
#'     \eqn{V(\mu) = \sigma^2} is independent of the mean and the distinction
#'     does not arise.}
#'   \item{temporal_var}{(Only present when the model formula contains an
#'     autocorrelation term such as \code{ar()}, \code{ma()},
#'     \code{arma()}, \code{cosy()}, \code{unstr()}, \code{sar()}, or
#'     \code{car()}.) Variance attributable to residual temporal or spatial
#'     dependence beyond the named channels, computed as the signed remainder
#'     \code{total_var - (param_var + env_var + residual_var)}.  This captures
#'     the AR/MA/ARMA innovation accumulation that a single-innovation
#'     \code{residual_var} cannot represent.  It is \strong{signed}: a
#'     negative value means autocorrelated errors \emph{reduce} predictive
#'     variance relative to the independent-innovation baseline, which is a
#'     real possibility for a covariance term.  Being a remainder, it also
#'     absorbs Monte Carlo error from the other channels, so read a
#'     \code{temporal_var} of the same order as \code{v_env_mcse} as noise.}
#'   \item{closure_gap}{(Only present for iid models.)
#'     \code{total_var - (param_var + env_var + residual_var)}: how far the
#'     named channels fall short of the measured posterior predictive variance.
#'     It is Monte Carlo error for a correctly specified iid model, of the order
#'     of \code{v_env_mcse}.  A systematically large value means a channel is
#'     mis-estimated.}
#'   \item{total_var}{Total predictive variance: the \strong{sampled variance
#'     of the environmentally inclusive posterior predictive draws}, which is
#'     the distribution the reported credible interval is drawn from.  The
#'     driver perturbation is applied \emph{before} the predictive draw, so
#'     \code{env_var} is a genuine share of this total.}
#' }
#'
#' Changed in 1.4.0.  \code{total_var} used to be \emph{defined} as the sum of
#' the other components, which made the shares reconcile to 100\% by
#' construction and made \code{env_var} an add-on to a predictive variance that
#' did not contain it.  It is now measured, and the sum is recovered as a
#' \emph{result}: the channels add up to \code{total_var} exactly for
#' autocorrelation models (where \code{temporal_var} is the remainder) and up to
#' Monte Carlo error for iid models (reported as \code{closure_gap}).
#' \code{param_var}, \code{env_var} and \code{residual_var} remain
#' non-negative; \code{temporal_var} and \code{closure_gap} are signed.
#'
#' @param predictions An \code{et_prediction} object from
#'   \code{\link{et_predict}}, or an \code{et_prediction_list} (grouped).
#' @param ... Unused.
#' @return A \code{data.frame} with columns
#'   \code{obs_id, param_var, env_var, residual_var, total_var}
#'   (plus \code{temporal_var} for autocorrelation models or
#'   \code{closure_gap} for iid ones, and a leading \code{group} column for
#'   grouped predictions).
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
