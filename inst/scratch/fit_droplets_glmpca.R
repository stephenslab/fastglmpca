command_args = commandArgs(trailingOnly = TRUE)
n_factor = as.integer(command_args[1])
n_iter = as.integer(command_args[2])
optimizer = as.character(command_args[3])
minibatch = as.character(command_args[4])

print(optimizer)
print(minibatch)

load("/project2/mstephens/pcarbo/git/fastTopics-experiments/data/droplet.RData")
#load("~/Documents/plash/scratch/droplet.RData")

if (optimizer == "fisher") {
  
  data <- as.matrix(counts)
  
}


set.seed(1)
if (optimizer == "fisher") {
  
  fit0 <- plash::init_glmpca(
    Y = data, K = n_factor, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
  )
  
} else {
  
  fit0 <- plash::init_glmpca(
    Y = counts, K = n_factor, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
  )
  
}


library(tictoc)

tic()
if(optimizer == "fisher") {
  
  fit <- glmpca::glmpca(
    Y = data,
    L = n_factor,
    fam = "poi",
    optimizer = optimizer,
    minibatch = minibatch,
    ctl = list(minIter = 1, maxIter = n_iter, verbose = TRUE, tol = .Machine$double.eps),
    init = list(factors = t(fit0$FF[-c(1,2),]), loadings = t(fit0$LL[-c(1,2),]))
  )
  
} else {
  
  fit <- glmpca::glmpca(
    Y = counts,
    L = n_factor,
    fam = "poi",
    optimizer = optimizer,
    minibatch = minibatch,
    ctl = list(minIter = 11, maxIter = n_iter, verbose = TRUE, tol = .Machine$double.eps),
    init = list(factors = t(fit0$FF[-c(1,2),]), loadings = t(fit0$LL[-c(1,2),]))
  )
  
}

toc()

readr::write_rds(
  fit, 
  glue::glue("droplets_glmpca_fit_{n_factor}_factors_{n_iter}_iter_{optimizer}_optimizer_minibatch_{minibatch}.rds")
)