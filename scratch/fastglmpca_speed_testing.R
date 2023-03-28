# here, want to design a test of computational speed for glmpca data

n_cells <- c(1000)
n_genes <- 2 * n_cells

time_vec <- c()
loglik_vec <- c()

for (i in 1:length(n_cells)) {
  
  set.seed(n_cells[i])
  sim_data <- plash::generate_glmpca_data(
    n = n_genes[i], p = n_cells[i], K = 8
  )
  
  sim_data$Y <- sim_data$Y[rowSums(sim_data$Y) > 0,]
  sim_data$Y <- sim_data$Y[, colSums(sim_data$Y) > 0]
  
  set.seed(1)
  
  fit0 <- plash::init_glmpca(
    Y = sim_data$Y, K = 8, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
  )
  tictoc::tic()
  fit <- plash::fit_glmpca(Y = sim_data$Y, fit0 = fit0, control = list(line_search = TRUE, num_iter = 5))
  stop <- tictoc::toc()
  time_elapsed <- stop$toc - stop$tic
  loglik <- tail(fit$progress$loglik, 1)
  
  time_vec <- c(time_vec, time_elapsed)
  loglik_vec <- c(loglik_vec, loglik)
  
}

res_df <- data.frame(
  time = time_vec, loglik = loglik_vec, ncells = n_cells, ngenes = n_genes
)

readr::write_csv(
  res_df, "sim_results_fastGLMPCA.csv"
)
