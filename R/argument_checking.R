# Verify that x is matrix with finite, numeric entries.
verify.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a numeric matrix, and",
               "all entries should be finite and non-missing")
  if (!(is.matrix(x) & is.numeric(x)) & !inherits(x, "sparseMatrix"))
    stop(msg)
  else if (any(is.infinite(x)) | any(is.na(x)))
    stop(msg)
  return(TRUE)
}

# Verify that x is a valid count matrix.
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
verify.count.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a non-negative,",
               "numeric matrix with at least 2 rows and 2 columns,",
               "all entries should be finite and non-missing, and",
               "no rows or columns should be entirely zeros")
  verify.matrix(x,arg.name)
  if (!(nrow(x) > 1 & ncol(x) > 1))
    stop(msg)
  else if (any(x < 0))
    stop(msg)
  else if (!(all(rowSums(x) > 0) & all(colSums(x) > 0)))
    stop(msg)
  return(TRUE)
}

# Return TRUE if x is a finite scalar with no missing entries.
is.scalar <- function (x)
  is.numeric(x) &
  length(x) == 1 &
  all(!is.na(x)) &
  all(is.finite(x))

# Verify that x is a valid GLM-PCA fit.
verify.fit <- function (x, arg.name = deparse(substitute(x))) {
  arg.name.LL <- paste0(arg.name,"$LL")
  arg.name.FF <- paste0(arg.name,"$FF")
  arg.name   <- sprintf("\"%s\"",arg.name)
  msg        <- paste("Input argument",arg.name,"should be a list containing",
                      "matrices \"LL\" and \"FF\"")
  if (!is.list(x))
    stop(msg)
  else if (!all(is.element(c("LL","FF"),names(x))))
    stop(msg)
  verify.matrix(x$LL,arg.name.LL)
  verify.matrix(x$FF,arg.name.FF)
  arg.name.LL <- sprintf("\"%s\"",arg.name.LL)
  arg.name.FF <- sprintf("\"%s\"",arg.name.FF)
  if (nrow(x$LL) != nrow(x$FF))
    stop(paste("Input matrices",arg.name.LL,"and",arg.name.FF,"should have",
               "the same number of rows"))
  return(TRUE)
}

