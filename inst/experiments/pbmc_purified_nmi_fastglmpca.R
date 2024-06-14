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

# while the total time is less than 10 hours
while (total_time < (10 * 60 * 60)) {
  
  d <- distances(fit0$V)
  dm <- distance_matrix(d)
  clust_tree <- fastcluster::hclust(dm, method="ward.D2")
  clusts <- cutree(clust_tree, k = 10)
  nmi <- NMI(celltype, clusts)
  nmi_vec <- c(nmi_vec, nmi)
  
  start_iter_time <- Sys.time()
  # update fit
  fit0 <- fit_glmpca_pois(
    Y = counts,
    fit0 = fit0,
    control = list(maxiter = 2)
  )
  end_iter_time <- Sys.time()
  total_time <- total_time + as.numeric(
    difftime(end_iter_time, start_iter_time, units = "secs")
  )
  
}

res_df <- data.frame(
  iter = 0:(length(nmi_vec) - 1),
  nmi = nmi_vec
)

readr::write_rds(res_df, "pbmc_purified_nmi_fastglmpca_10hr_res_k10_ward.rds")
readr::write_rds(fit0, "pbmc_purified_k10_10hr_iter.rds")

