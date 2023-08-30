command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])
n_cores = as.character(command_args[3])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/pbmc_68k.RData")

library(RhpcBLASctl)
omp_set_num_threads(n_cores)
blas_set_num_threads(1)

set.seed(1) 
fit0 <- fastglmpca::init_glmpca_pois(
  Y = Matrix::t(counts), 
  K = n_factor, 
  fit_col_size_factor = TRUE,
  fit_row_intercept = TRUE
)

fit <- fastglmpca::fit_glmpca_pois(
  Y = Matrix::t(counts),
  fit0 = fit0,
  min_iter = n_iter - 1,
  max_iter = n_iter
)

readr::write_rds(
  fit,
  glue::glue("pbmc_fastglmpca_fit_{n_factor}_factors_{n_iter}_iter_{n_cores}_cores.rds")
)