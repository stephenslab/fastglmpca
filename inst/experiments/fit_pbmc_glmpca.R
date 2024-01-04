command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/pbmc_68k.RData")

set.seed(1)

fit0 <- fastglmpca::init_glmpca_pois(
  Y = Matrix::t(counts),
  K = n_factor
)

library(tictoc)

fit <- glmpca::glmpca(
  Y = Matrix::t(counts),
  L = n_factor,
  fam = "poi",
  optimizer = optimizer,
  minibatch = minibatch,
  ctl = list(
    minIter = n_iter - 1,
    maxIter = n_iter,
    verbose = TRUE,
    tol = .Machine$double.eps,
    lr = 1e-4
  ),
  init = list(
    factors = fit0$V %*% diag(sqrt(fit0$d)),
    loadings = fit0$U %*% diag(sqrt(fit0$d))
  )
)

readr::write_rds(
  fit,
  glue::glue("pbmc_glmpca_fit_{n_factor}_factors_{n_iter}_iter_avagrad_optimizer_minibatch_stochastic_dec_23.rds")
)
