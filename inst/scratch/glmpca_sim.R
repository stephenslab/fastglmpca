library(fastTopics)
library(Matrix)

get_lik_nmf <- function(Y, ft_fit) {
  
  Lambda <- tcrossprod(ft_fit$L, ft_fit$F)
  lik <- dpois(drop(Y), drop(Lambda), log = FALSE)
  return(mean(lik))
  
}

simulate_scRNA_data <- function(n_cells, n_genes, K, method = c("glmpca", "nmf")) {
  
  if (method == "nmf") {
    
    sim <- fastTopics::simulate_multinom_gene_data(n = n_cells, m = n_genes, k = K)
    return(sim$X)
    
  } else if (method == "glmpca") {
    
    sim <- plash:::generate_plash_data(n = n_genes, p = n_genes, k = K)
    return(t(sim$Y))
    
  }
  
}

# Now, try and fit my algo
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"
)

doParallel::registerDoParallel(cl = my.cluster)

n_sims_per_iter <- 10
n_cells_vec <- c()
n_genes_vec <- c()
K_vec <- c()
glmpca_ll_vec <- c()
nmf_ll_vec <- c()

for (n_cells in c(250, 500, 1000, 2000, 3500)) {
  
  n_genes <- 2 * n_cells
  for (K in c(2, 5, 10, 20)) {
    
    avg_glmpca_ll <- 0
    avg_nfm_ll <- 0
    
    for (i in 1:n_sims_per_iter) {
      
      nmf_data <- simulate_scRNA_data(
        n_cells = n_cells, n_genes = n_genes, K = K, method = "glmpca"
      )
      
      nmf_fit <- fastTopics::fit_poisson_nmf(as(nmf_data, "sparseMatrix"), k = K)
      glmpca_fit <- plash::plash_omni(
        Y = t(nmf_data),
        K = K, 
        offset = TRUE,
        intercept = FALSE,
        update_c = FALSE,
        init_cc = rep(0, n_genes),
        parallel = TRUE,
        tol = 5e-5,
        min_iter = 3
      )
      
      nmf_lik <- get_lik_nmf(nmf_data, nmf_fit)
      glmpca_lik <- glmpca_fit$lik / (n_genes * n_cells)
      
      avg_glmpca_ll <- avg_glmpca_ll + (1 / n_sims_per_iter) * glmpca_lik
      avg_nfm_ll <- avg_nfm_ll + (1 / n_sims_per_iter) * nmf_lik
      
    }
    
    n_cells_vec <- c(n_cells_vec, n_cells)
    n_genes_vec <- c(n_genes_vec, n_genes)
    K_vec <- c(K_vec, K)
    glmpca_ll_vec <- c(glmpca_ll_vec, avg_glmpca_ll)
    nmf_ll_vec <- c(nmf_ll_vec, avg_nfm_ll)
    
  }
  
}

sim_df <- data.frame(
  n_genes = n_genes_vec,
  n_cells = n_cells_vec,
  K = K_vec,
  glmpca_loglik = glmpca_ll_vec,
  nmf_loglik = nmf_ll_vec
)

readr::write_csv(sim_df, "./glmpca_data_sim_res.csv")