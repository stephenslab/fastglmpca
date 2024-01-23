set.seed(1)

ncells_vec <- c(100, seq(200, 1000), 1500, 2500, 5000)
nsims_per_exp <- 100

fastglmpca_time_vec <- c()
glmpca_time_vec <- c()
scGBM_time_vec <- c()

for (ncells in ncells_vec) {
  
  fastglmpca_avg_time <- 0
  glmpca_avg_time <- 0
  scGBM_avg_time <- 0
  
  for (i in 1:nsims_per_exp) {
    
    print(
      glue::glue("Printing results for simulation {i} with {ncells} cells")
    )
    
    sim_data <- fastglmpca::generate_glmpca_data_pois(
      n = 2000, p = ncells, K = 5
    )
    
    counts <- as(sim_data$Y, "sparseMatrix")
    counts <- counts[(Matrix::rowSums(counts) > 0), ]
    counts <- counts[, (Matrix::colSums(counts) > 0)]
    
    fit0 <- fastglmpca::init_glmpca_pois(
      Y = counts,
      K = 5
    )
    
    fastglmpca_fit <- fastglmpca::fit_glmpca_pois(
      Y = counts,
      fit0 = fit0,
      control = list(maxiter = 10)
    )
    fastglmpca_finish <- sum(fastglmpca_fit$progress$time)
    
    if (i == 2 && ncells == 200) {
      
      print(0)
      
    }
    
    glmpca_fit <- glmpca::glmpca(
      Y = counts,
      L = 5,
      fam = "poi",
      optimizer = "avagrad",
      minibatch = "stochastic",
      ctl = list(
        maxIter = 10,
        minIter = 9,
        lr = 1e-4,
        batch_size = ceiling(ncells / 10)
      )
    )
    glmpca_finish <- sum(glmpca_fit$time)
    
    scGBM_fit <- scGBM::gbm.sc(
      Y = as.matrix(counts),
      M = 5,
      max.iter = 10,
      min.iter = 10,
      return.W = FALSE,
      time.by.iter = TRUE
    )
    scGBM_finish <- sum(scGBM_fit$time)
    
    fastglmpca_avg_time <- fastglmpca_avg_time + (1 / nsims_per_exp) * fastglmpca_finish
    glmpca_avg_time <- glmpca_avg_time + (1 / nsims_per_exp) * glmpca_finish
    scGBM_avg_time <- scGBM_avg_time + (1 / nsims_per_exp) * scGBM_finish
    
  }
  
  fastglmpca_time_vec <- c(fastglmpca_time_vec, fastglmpca_avg_time)
  glmpca_time_vec <- c(glmpca_time_vec, glmpca_avg_time)
  scGBM_time_vec <- c(scGBM_time_vec, scGBM_avg_time)
  
}

time_n_cells_df <- data.frame(
  ncells = ncells_vec,
  fastglmpca_time = fastglmpca_time_vec,
  glmpca_time = glmpca_time_vec,
  scGBM_time = scGBM_time_vec
)

readr::write_rds(time_n_cells_df, "num_cells_scaling_df.rds")

# added 15 extra cell to account for cells that will be excluded since
# they have no counts > 0
ncells <- 1515
ngenes <- 2000

K_vec <- c(seq(2, 10), 15, 20, 25, 50)

nsims_per_exp <- 50

fastglmpca_time_vec <- c()
glmpca_time_vec <- c()
scGBM_time_vec <- c()

for (K in K_vec) {
  
  fastglmpca_avg_time <- 0
  glmpca_avg_time <- 0
  scGBM_avg_time <- 0
  
  for (i in 1:nsims_per_exp) {
    
    sim_data <- fastglmpca::generate_glmpca_data_pois(
      n = ncells, p = ngenes, K = 10
    )
    
    counts <- as(sim_data$Y, "sparseMatrix")
    counts <- counts[(Matrix::rowSums(counts) > 0), ]
    counts <- counts[, (Matrix::colSums(counts) > 0)]
    
    fit0 <- fastglmpca::init_glmpca_pois(
      Y = counts,
      K = K
    )
    
    fastglmpca_fit <- fastglmpca::fit_glmpca_pois(
      Y = counts,
      fit0 = fit0,
      control = list(maxiter = 10)
    )
    fastglmpca_finish <- sum(fastglmpca_fit$progress$time)
    
    glmpca_fit <- glmpca::glmpca(
      Y = counts,
      L = K,
      fam = "poi",
      optimizer = "avagrad",
      minibatch = "stochastic",
      ctl = list(
        maxIter = 10,
        minIter = 9,
        lr = 1e-4,
        batch_size = ceiling(ncells / 10)
      )
    )
    glmpca_finish <- sum(glmpca_fit$time)
    
    scGBM_fit <- scGBM::gbm.sc(
      Y = as.matrix(counts),
      M = K,
      max.iter = 10,
      min.iter = 10,
      return.W = FALSE,
      time.by.iter = TRUE
    )
    scGBM_finish <- sum(scGBM_fit$time)
    
    fastglmpca_avg_time <- fastglmpca_avg_time + (1 / nsims_per_exp) * fastglmpca_finish
    glmpca_avg_time <- glmpca_avg_time + (1 / nsims_per_exp) * glmpca_finish
    scGBM_avg_time <- scGBM_avg_time + (1 / nsims_per_exp) * scGBM_finish
    
  }
  
  fastglmpca_time_vec <- c(fastglmpca_time_vec, fastglmpca_avg_time)
  glmpca_time_vec <- c(glmpca_time_vec, glmpca_avg_time)
  scGBM_time_vec <- c(scGBM_time_vec, scGBM_avg_time)
  
}

time_k_cells_df <- data.frame(
  K = K_vec,
  fastglmpca_time = fastglmpca_time_vec,
  glmpca_time = glmpca_time_vec,
  scGBM_time = scGBM_time_vec
)

readr::write_rds(time_k_cells_df, "num_factors_scaling_df.rds")

