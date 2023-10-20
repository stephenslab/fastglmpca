#' @title Set up Multithreading for fastglmpca
#' 
#' @description Initialize RcppParallel multithreading using a
#'   pre-specified number of threads, or using the default number of
#'   threads when \code{n} is not specified or is NA.
#'
#' @param n The requested number of threads.
#'
#' @return The number of threads to be used.
#' 
#' @importFrom RcppParallel setThreadOptions
#' @importFrom RcppParallel defaultNumThreads
#'
#' @export
#' 
set_fastglmpca_threads <- function (n) {
  if (missing(n)) {
    setThreadOptions()
    n <- defaultNumThreads()
  } else if (is.na(n)) {
    setThreadOptions()
    n <- defaultNumThreads()
  }
    setThreadOptions(numThreads = n)
  message(sprintf("Using %d RcppParallel threads.",n))
  return(n)
}
