get_cluster_labels <- function(V) {
  
  d <- distances::distances(V)
  dm <- distances::distance_matrix(d)
  clust_tree <- fastcluster::hclust(dm, method="ward.D2")
  clusts <- cutree(clust_tree, k = 10)
  names(clusts) <- rownames(V)
  return(clusts)
  
}

fit_fastglmpca28 <- readr::read_rds("pbmc_purified_k10_10hr_iter.rds")
fit_fastglmpca1 <- readr::read_rds("pbmc_purified_k10_10hr_iter_1core.rds")
fit_scGBM <- readr::read_rds("pbmc_purified_scGBM_k10_10hr_iter.rds")
fit_glmpca <- readr::read_rds("pbmc_purified_glmpca_k10_10hr_iter.rds")

rownames(fit_glmpca$V) <- rownames(fit_fastglmpca28$V)

clusters_fastglmpca28 <- get_cluster_labels(fit_fastglmpca28$V)
clusters_fastglmpca1 <- get_cluster_labels(fit_fastglmpca1$V)
clusters_scGBM <- get_cluster_labels(fit_scGBM$V)
clusters_glmpca <- get_cluster_labels(fit_glmpca$V)

save(
  clusters_fastglmpca28,
  clusters_fastglmpca1,
  clusters_scGBM,
  clusters_glmpca,
  file = "pbmc_purified_clusters.Rdata"
)

