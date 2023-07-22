command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/pbmc_68k.RData")

data <- as.matrix(Matrix::t(counts))

set.seed(1)

library(tictoc)

tic()
fit <- scGBM::gbm.sc(
  Y = data,
  M = n_factor,
  max.iter = 1e8,
  tol = .Machine$double.eps,
  time.by.iter = TRUE,
  infer.beta = FALSE,
  return.W = FALSE,
  min.iter = 1e8 - 1
)
toc()

readr::write_rds(
  fit,
  glue::glue("pbmc_scGBM_fit_{n_factor}_factors_no_beta_infer_10_hrs.rds")
)