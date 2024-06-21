load("/project2/mstephens/pcarbo/git/single-cell-topics/data/pbmc_purified.RData")
library(Matrix)
set.seed(1)
counts <- counts[, Matrix::colSums(counts) > 0]
counts <- as.matrix(Matrix::t(counts))
gc()

celltype <- samples$celltype

fit <- scGBM::gbm.sc(
  Y = counts,
  M = 10,
  max.iter = 1e6,
  return.W = FALSE,
  tol = .Machine$double.eps,
  time.by.iter = TRUE,
  cluster = TRUE,
  celltype = celltype,
  max_seconds_run = 10 * 60 * 60
)

readr::write_rds(fit, "scGBM_10hours_pbmc_purified_nmi_ari.rds")
