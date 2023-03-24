# here, want to design a test of computational speed for glmpca data

n_cells <- c(1000, 2500, 5000, 10000, 25000)
n_genes <- 2 * n_cells

time_vec <- c()
loglik_vec <- c()

for (i in 1:length(n_cells)) {
  
  set.seed(n_cells)
  sim_data <- plash::generate_glmpca_data(
    n = n_genes, p = n_cells, K = 10
  )
  
  sim_data$Y <- sim_data$Y[rowSums(data$Y) > 0,]
  sim_data$Y <- sim_data$Y[, colSums(data$Y) > 0]
  
  set.seed(1)
  
  tictoc::tic()
  fit0 <- plash::init_glmpca(
    Y = sim_data$Y, K = 10, fit_col_size_factor = TRUE, fit_row_intercept = TRUE
  )
  stop <- tictoc::toc()
  time_elapsed <- stop$toc - stop$tic
  loglik <- tail(fit0$progress$loglik, 1)
  
  time_vec <- c(time_vec, time_elapsed)
  loglik_vec <- c(loglik_vec, loglik)
  
}

res_df <- data.frame(
  time = time_vec, loglik = loglik_vec, ncells = n_cells, ngenes = n_genes
)

readr::write_csv(
  res_df, "sim_results_fastGLMPCA.csv"
)
