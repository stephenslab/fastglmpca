# Used to check that x is a vector in which x[i+1] >= x[i] for all i.
expect_nondecreasing <- function (x, tol=1e-7)
  expect_equal(diff(x) >= -tol,rep(TRUE,length(x) - 1))
