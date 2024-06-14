load("/project2/mstephens/pcarbo/git/single-cell-topics/data/pbmc_purified.RData")
#load("~/Documents/data/fastglmpca/raw_data/pbmc_purified.RData")
library(Matrix)
set.seed(1)
counts <- counts[, Matrix::colSums(counts) > 0]
counts <- as.matrix(Matrix::t(counts))
gc()

library(aricode)
library(fastglmpca)
library(igraph)
library(cccd)
library(distances)
library(scGBM)
#library(Rclusterpp)
celltype <- samples$celltype  #[samp_idx]
rm(samples, genes)
gc()

set.seed(1)
start_iter_time <- Sys.time()
fit0 <- gbm.init(Y = counts, M = 10)
end_iter_time <- Sys.time()
total_time <- as.numeric(difftime(end_iter_time, start_iter_time, units="secs"))

nmi_vec <- c()

#total_time <- 0

iter <- 0

# while the total time is less than 10 hours
while (total_time < (10 * 60 * 60)) {
  
  print(iter)
  iter <- iter + 1
  
  start_clust_time <- Sys.time()
  
  d <- distances(fit0$V)
  dm <- distance_matrix(d)
  clust_tree <- fastcluster::hclust(dm, method="ward.D2")
  #clust_tree <- hclust(dm, method="ward.D2")
  clusts <- cutree(clust_tree, k = 10)
  nmi <- NMI(celltype, clusts)
  nmi_vec <- c(nmi_vec, nmi)
  
  end_clust_time <- Sys.time()
  clust_time <- as.numeric(difftime(end_clust_time, start_clust_time, units="secs"))
  print("cluster time:")
  print(clust_time)
  
  start_iter_time <- Sys.time()
  
  if (iter == 1) {
    
    # update fit with glmpca
    fit0 <- gbm.sc(
      Y = counts,
      M = 10,
      init_U = fit0$U,
      init_V = fit0$V,
      init_d = fit0$D,
      max.iter = 1,
      tol = .Machine$double.eps,
      time.by.iter = FALSE,
      infer.beta = FALSE,
      return.W = FALSE
    )
    
  } else {
    
    # update fit with glmpca
    fit0 <- gbm.sc(
      Y = counts,
      M = 10,
      init_U = fit0$U,
      init_V = fit0$V,
      init_d = fit0$D,
      max.iter = 1,
      tol = .Machine$double.eps,
      time.by.iter = FALSE,
      infer.beta = FALSE,
      return.W = FALSE,
      lr = fit0$lr,
      LL = fit0$LL,
      current_iter = iter
    )
    
  }
  

  
  end_iter_time <- Sys.time()
  total_time <- total_time + as.numeric(
    difftime(end_iter_time, start_iter_time, units = "secs")
  )
  
  print("total time:")
  print(total_time)
  
}

res_df <- data.frame(
  iter = 0:length(nmi_vec),
  nmi = c(0, nmi_vec)
)

readr::write_rds(res_df, "pbmc_purified_nmi_scGBM_10hr_res_k10_ward.rds")
readr::write_rds(fit0, "pbmc_purified_scGBM_k10_10hr_iter.rds")
