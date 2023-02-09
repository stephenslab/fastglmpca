#' @rdname fitted.glmpcafit
#'
#' @title Get Fitted Values for GLM-PCA Model Fit
#'
#' @description \code{fitted} method for the
#'   \dQuote{glmpca_fit} class.
#'
#' @param object An object of class \dQuote{glmpca_fit},
#'   typically the result of calling \code{\link{fit_glmpca}}.
#'
#' @return An n x p matrix of fitted means. Calculated as
#'   \deqn{exp(L'F)} using the \code{fit} object.
#' 
#' @method fitted glmpca_fit
#'
#' @export
#' 
fitted.glmpca_fit <- function (object, ...) {
  verify.fit(object)
  exp(crossprod(object$LL, object$FF))
}

#' @rdname summary.glmpca_fit
#'
#' @title Summarize GLM-PCA Model Fit
#'
#' @description \code{summary} method for the
#'   \dQuote{glmpcan_fit} class.
#'
#' @param object An object of class \dQuote{glmpca_fit},
#'   typically the result of calling \code{\link{fit_glmpca}}.
#'
#' @param x An object of class \dQuote{summary.glmpca_fit},
#'   usually the result of a call to \code{summary.glmpca_fit}.
#'
#' @param \dots Additional arguments passed to the generic
#'   \code{summary} or \code{print.summary} method.
#'
#' @return A vector of statistics summarizing the model fit.
#' 
#' @method summary glmpca_fit
#'
#' @export
#' 
summary.glmpca_fit <- function (object, ...) {
  verify.fit(object)
  numiter <- nrow(object$progress)
  out <- c(n       = ncol(object$LL),
           p       = ncol(object$FF),
           K       = nrow(object$LL),
           numiter = numiter,
           loglik  = object$progress[numiter,"loglik"])
  class(out) <- c("summary.glmpca_fit","list")
  return(out)
}

#' @rdname summary.glmpca_fit
#'
#' @method print summary.glmpca_fit
#'
#' @export
#'
print.summary.glmpca_fit <- function (x, ...) {
  cat("Model overview:\n")
  cat(sprintf("  Number of data rows, n: %d\n",x["n"]))
  cat(sprintf("  Number of data cols, p: %d\n",x["p"]))
  cat(sprintf("  Rank: %d\n",x["K"]))
  cat(sprintf("Evaluation of model fit (%d updates performed):\n",
              x["numiter"]))
  cat(sprintf("  log-likelihood: %+0.12e\n",x["loglik"]))
  return(invisible(x))
}