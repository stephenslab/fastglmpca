command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])
optimizer = as.character(command_args[4])
minibatch = as.character(command_args[5])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")

data <- as.matrix(counts)

set.seed(1)
fit0 <- plash::init_glmpca(
  Y = data, K = n_factor, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
)

tic()
fit <- glmpca::glmpca(
  Y = data,
  L = n_factors,
  fam = "poi",
  optimizer = optimizer,
  ctl = list(minIter = 1, maxIter = n_iter, verbose = TRUE, tol = .Machine$double.eps),
  init = list(factors = t(fit0$FF[-c(1,2),]), loadings = t(fit0$LL[-c(1,2),]))
)
toc()

readr::write_rds(
  fit, 
  glue::glue("droplets_glmpca_fit_{n_factor}_factors_{n_iter}_iter_{optimizer}_optimizer_minibatch_{minibatch}.rds")
)