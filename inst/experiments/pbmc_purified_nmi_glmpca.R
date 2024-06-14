load("/project2/mstephens/pcarbo/git/single-cell-topics/data/pbmc_purified.RData")
library(Matrix)
set.seed(1)
counts <- counts[, Matrix::colSums(counts) > 0]
counts <- Matrix::t(counts)
gc()

library(aricode)
library(fastglmpca)
library(igraph)
library(cccd)
library(distances)
#library(Rclusterpp)
celltype <- samples$celltype  #[samp_idx]
rm(samples, genes)
gc()

set.seed(1)
fit0 <- init_glmpca_pois(
  Y = counts,
  K = 10
)

nmi_vec <- c()

total_time <- 0

iter <- 0

# while the total time is less than 10 hours
while (total_time < (10 * 60 * 60)) {
  
  print(iter)
  iter <- iter + 1
  
  d <- distances(fit0$V)
  dm <- distance_matrix(d)
  clust_tree <- fastcluster::hclust(dm, method="ward.D2")
  #clust_tree <- hclust(dm, method="ward.D2")
  clusts <- cutree(clust_tree, k = 10)
  nmi <- NMI(celltype, clusts)
  nmi_vec <- c(nmi_vec, nmi)
  
  start_iter_time <- Sys.time()
  # update fit with glmpca
  if (iter > 1) {
    
    fit_glmpca <- glmpca::glmpca(
      Y = counts,
      L = 10,
      fam = "poi",
      calc_likelihood = FALSE,
      optimizer = "avagrad",
      minibatch = "stochastic",
      m_u=fit_glmpca$m_u,
      m_v=fit_glmpca$m_v,
      v_u=fit_glmpca$v_u,
      v_v=fit_glmpca$v_v,
      ctl = list(
        minIter = 0,
        maxIter = 1,
        verbose = TRUE,
        tol = .Machine$double.eps,
        lr = 5e-5
      ),
      init = list(
        factors = fit0$V %*% diag(sqrt(fit0$d)),
        loadings = fit0$U %*% diag(sqrt(fit0$d))
      )
    )
    
  }
  
  else {
    
    fit_glmpca <- glmpca::glmpca(
      Y = counts,
      L = 10,
      fam = "poi",
      calc_likelihood = FALSE,
      optimizer = "avagrad",
      minibatch = "stochastic",
      ctl = list(
        minIter = 0,
        maxIter = 1,
        verbose = TRUE,
        tol = .Machine$double.eps,
        lr = 5e-5
      ),
      init = list(
        factors = fit0$V %*% diag(sqrt(fit0$d)),
        loadings = fit0$U %*% diag(sqrt(fit0$d))
      )
    )
    
  }
  
  tmp_mod <- fastglmpca:::orthonormalize(
    as.matrix(fit_glmpca$loadings), as.matrix(fit_glmpca$factors)
  )
  
  end_iter_time <- Sys.time()
  total_time <- total_time + as.numeric(
    difftime(end_iter_time, start_iter_time, units = "secs")
  )
  
  fit0$U <- tmp_mod$U
  fit0$V <- tmp_mod$V
  fit0$d <- tmp_mod$d
  
}

res_df <- data.frame(
  iter = 0:(length(nmi_vec) - 1),
  nmi = nmi_vec
)

readr::write_rds(res_df, "pbmc_purified_nmi_glmpca_10hr_res_k10_ward.rds")
readr::write_rds(fit0, "pbmc_purified_glmpca_k10_10hr_iter.rds")
