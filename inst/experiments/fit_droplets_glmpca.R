command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])

library(Matrix)
library(RhpcBLASctl)
blas_set_num_threads(1)
omp_set_num_threads(1)
load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")

set.seed(1)
fit0 <- fastglmpca::init_glmpca_pois(
  Y = Matrix::t(counts), K = n_factor, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
)

library(tictoc)

tic()
fit <- glmpca::glmpca(
  Y = Matrix::t(counts),
  L = n_factor,
  fam = "poi",
  optimizer = "avagrad",
  minibatch = "stochastic",
  ctl = list(minIter = 1e8 - 1, maxIter = 1e8, verbose = TRUE, tol = .Machine$double.eps, lr = 4e-11),
  init = list(factors = fit0$V[, -c(1,2)], loadings = fit0$U[, -c(1,2)])
)

toc()

readr::write_rds(
  fit,
  glue::glue(
    "droplets_glmpca_fit_{n_factor}_factors_avagrad_optimizer_minibatch_stochastic_10_hrs.rds"
  )
)