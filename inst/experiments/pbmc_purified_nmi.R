#load("/project2/mstephens/pcarbo/git/single-cell-topics/data/pbmc_purified.RData")
load("~/Documents/data/fastglmpca/raw_data/pbmc_purified.RData")

counts <- counts[, Matrix::colSums(counts) > 0]
counts <- Matrix::t(counts)
gc()

celltype <- samples$celltype
rm(samples, genes)
gc()

set.seed(1)
fit0 <- init_glmpca_pois(
  Y = counts,
  K = 25
)

nmi_vec <- c()

for (iter in 1:51) {
  
  print(iter)
  
  # first, get clustering based on previous model
  kmm <- kmeans(
    x = fit0$V,
    centers = 7,
    nstart = 5,
    iter.max = 25
  )
  nmi <- NMI(celltype, kmm$cluster)
  nmi_vec <- c(nmi_vec, nmi)
  
  # update fit
  fit0 <- fit_glmpca_pois(
    Y = counts, 
    fit0 = fit0,
    control = list(maxiter = 2)
  )
  
}

res_df <- data.frame(
  iter = 0:50,
  nmi = nmi_vec
)

readr::write_rds(res_df, "pbmc_purified_nmi_res.rds")
readr::write_rds(fit0, "pbmc_purified_k25_50iter.rds")


