# Return true if x is a compressed, sparse, column-oriented numeric
# matrix.
is.sparse.matrix <- function (x)
  inherits(x,"dsparseMatrix") && is.numeric(x@x)

# Verify that x is matrix with finite, numeric entries.
verify.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a numeric matrix, and",
               "all entries should be finite and non-missing")
  if (!(is.matrix(x) & is.numeric(x)) & !is.sparse.matrix(x))
    stop(msg)
  else if (any(is.infinite(x)) | anyNA(x))
    stop(msg)
  return(TRUE)
}

# Verify that x is a valid count matrix.
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' 
verify.count.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a non-negative,",
               "numeric matrix with at least 2 rows and 2 columns,",
               "and all entries should be finite and non-missing")
  verify.matrix(x,arg.name)
  if (!(nrow(x) > 1 & ncol(x) > 1))
    stop(msg)
  else if (any(x < 0))
    stop(msg)
  
  if (any(Matrix::rowSums(x) == 0)) 
    stop(
      paste(
        "Input argument",arg.name,"cannot have any rows with all 0 entries"
      )
    )
  if (any(Matrix::colSums(x) == 0)) 
    stop(
      paste(
        "Input argument",arg.name,"cannot have any columns with all 0 entries"
      )
    )
  
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
  arg.name.U <- paste0(arg.name,"$U")
  arg.name.V <- paste0(arg.name,"$V")
  arg.name   <- sprintf("\"%s\"",arg.name)
  msg        <- paste("Input argument",arg.name,"should be a list containing",
                      "matrices \"U\" and \"V\"")
  if (!is.list(x))
    stop(msg)
  else if (!all(is.element(c("U","V"),names(x))))
    stop(msg)
  verify.matrix(x$U,arg.name.U)
  verify.matrix(x$V,arg.name.V)
  arg.name.U <- sprintf("\"%s\"",arg.name.U)
  arg.name.V <- sprintf("\"%s\"",arg.name.V)
  if (ncol(x$U) != ncol(x$V))
    stop(paste("Input matrices",arg.name.U,"and",arg.name.V,"should have",
               "the same number of columns"))
  return(TRUE)
}

