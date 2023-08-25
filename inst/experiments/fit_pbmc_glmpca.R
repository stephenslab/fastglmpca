command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])
optimizer = as.character(command_args[3])
minibatch = as.character(command_args[4])

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/pbmc_68k.RData")

if (optimizer == "fisher") {
  
  data <- as.matrix(counts)
  
}

set.seed(1)
if (optimizer == "fisher") {
  
  fit0 <- fastglmpca::init_glmpca(
    Y = t(data), 
    K = n_factor, 
    fit_col_size_factor = TRUE,
    fit_row_intercept = TRUE
  )
  
} else if (optimizer == "avagrad"){
  
  fit0 <- fastglmpca::init_glmpca_pois(
    Y = Matrix::t(counts), 
    K = n_factor, 
    fit_col_size_factor = TRUE,
    fit_row_intercept = TRUE
  )
  
}


library(tictoc)

if(optimizer == "fisher") {
  
  fit <- glmpca::glmpca(
    Y = t(data),
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
    init = list(factors = fit0$V, loadings = fit0$U)
  )
  
} else {
  
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
    init = list(factors = fit0$V, loadings = fit0$U)
  )
  
}

readr::write_rds(
  fit,
  glue::glue("pbmc_glmpca_fit_{n_factor}_factors_{n_iter}_iter_{optimizer}_optimizer_minibatch_{minibatch}.rds")
)