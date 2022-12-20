# Again, I want to have a good understanding of when the NMF model
# will perform better than the GLMPCA model

# Once I know that, these becomes candidates for when my model will perform better 
# than GLMPCA

# I'm not too worried about my model being correct, as I've verified it repeatedly

get_lik_glmpca <- function(Y, glmpca_fit) {
  
  Lambda <- exp(crossprod(glmpca_fit$LL, glmpca_fit$FF)) - outer(glmpca_fit$cc, glmpca_fit$size)
  Lambda <- pmax(Lambda, .Machine$double.eps)
  lik <- dpois(drop(Y), drop(Lambda), log = TRUE)
  return(mean(lik))
  
}

get_lik_nmf <- function(Y, ft_fit) {
  
  Lambda <- tcrossprod(ft_fit$L, ft_fit$F)
  lik <- dpois(drop(Y), drop(Lambda), log = TRUE)
  return(mean(lik))
  
}

simulate_scRNA_data <- function(n_cells, n_genes, K, method = c("glmpca", "nmf")) {
  
  if (method == "nmf") {
    
    sim <- plash::simulate_poisson_gene_data(n = n_cells, m = n_genes, k = K, s =3, sparse = TRUE)
    return(sim)
    
  } else if (method == "glmpca") {
    
    # TODO
    
  }
  
}

# Now, try and fit my algo
n.cores <- parallel::detectCores()
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

n_sims_per_iter <- 5
n_cells_vec <- c()
n_genes_vec <- c()
K_vec <- c()
glmpca_ll_vec <- c()
nmf_ll_vec <- c()

for (n_cells in c(1000)) {
  
  n_genes <- ceiling(1.75 * n_cells)
  for (K in c(2, 5, 10, 20)) {
    
    avg_glmpca_ll <- 0
    avg_nfm_ll <- 0
    
    for (i in 1:n_sims_per_iter) {
      
      nmf_data <- simulate_scRNA_data(
        n_cells = n_cells, n_genes = n_genes, K = K, method = "nmf"
      )
      
      nmf_fit <- fastTopics::fit_poisson_nmf(nmf_data$X_train, k = K)
      glmpca_fit <- plash::plash_omni(
        Y = t(as.matrix(nmf_data$X_train)),
        K = K, 
        offset = FALSE,
        intercept = FALSE,
        update_c = FALSE,
        init_cc = rep(0, n_genes),
        parallel = TRUE,
        tol = 5e-5,
        min_iter = 3
      )
      
      nmf_lik <- get_lik_nmf(as.matrix(nmf_data$X_test), nmf_fit)
      glmpca_lik <- get_lik_glmpca(t(as.matrix(nmf_data$X_test)), glmpca_fit)
      
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

readr::write_csv(sim_df, "./nmf_data_sim_res.csv")
