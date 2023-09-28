#' @rdname fitted.glmpca_pois_fit
#'
#' @title Get Fitted Values for GLM-PCA Model Fit
#'
#' @description \code{fitted} method for the
#'   \dQuote{glmpca_pois_fit} class.
#'
#' @param object An object of class \dQuote{glmpca_fit},
#'   typically the result of calling \code{\link{fit_glmpca_pois}}.
#'   
#' @param \dots Additional arguments passed to the generic
#'   \code{fitted} method.
#'
#' @return An n x p matrix of fitted means. Calculated as
#'   \deqn{exp(UDV')} using the \code{fit} object.
#' 
#' @method fitted glmpca_pois_fit
#'
#' @export
#' 
fitted.glmpca_pois_fit <- function (object, ...) {
  verify.fit(object)
  return(exp(tcrossprod(object$U %*% object$D,object$V)))
}

#' @rdname summary.glmpca_pois_fit
#'
#' @title Summarize GLM-PCA Model Fit
#'
#' @description \code{summary} method for objects of class
#'   \dQuote{glmpcan_fit_pois}.
#'
#' @param object An object of class \dQuote{glmpca_fit},
#'   typically the result of calling \code{\link{fit_glmpca_pois}}.
#'
#' @param x An object of class \dQuote{summary.glmpca_fit},
#'   usually the result of a call to \code{summary.glmpca_fit}.
#'
#' @param \dots Additional arguments passed to the generic
#'   \code{summary} or \code{print.summary} method.
#'
#' @return \code{summary} returns a vector of basic statistics
#'   summarizing the model fit.
#' 
#' @method summary glmpca_pois_fit
#'
#' @export
#' 
summary.glmpca_pois_fit <- function (object, ...) {
  numiter <- length(object$progress$iter)
  out <- c(n       = nrow(object$U),
           m       = nrow(object$V),
           nx      = ifelse(length(object$X) > 0,ncol(object$X),0),
           nz      = ifelse(length(object$Z) > 0,ncol(object$Z),0),
           K       = ncol(object$U),
           numiter = max(object$progress$iter),
           loglik  = object$loglik)
  class(out) <- c("summary.glmpca_pois_fit","list")
  return(out)
}

#' @rdname summary.glmpca_pois_fit
#'
#' @method print summary.glmpca_pois_fit
#'
#' @export
#'
print.summary.glmpca_pois_fit <- function (x, ...) {
  cat(sprintf("GLM-PCA model fit to %d x %d count matrix:\n",x["n"],x["m"]))
  cat(sprintf("rank (K): %d\n",x["K"]))
  cat(sprintf("number of row covariates: %d\n",x["nx"]))
  cat(sprintf("number of column covariates: %d\n",x["nz"]))
  cat(sprintf("updates performed: %d\n",x["numiter"]))
  cat(sprintf("log-likelihood: %+0.8e\n",x["loglik"]))
  return(invisible(x))
}
