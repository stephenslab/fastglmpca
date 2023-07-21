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
fit <- fastglmpca::fit_glmpca_pois(
  Y = Matrix::t(counts),
  fit0 = fit0,
  control = list(line_search = TRUE, num_iter = 3, alpha = .01),
  max_iter = 1e8,
  min_iter = 1e8 - 1,
  tol = .Machine$double.eps
)
toc()

readr::write_rds(
  fit,
  glue::glue("droplets_fastglmpca_fit_1_core_{n_factor}_factors_10_hrs.rds")
)