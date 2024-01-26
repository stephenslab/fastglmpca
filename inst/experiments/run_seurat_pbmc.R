load("~/Documents/data/fastglmpca/raw_data/pbmc_68k.RData")

library(Seurat)
set.seed(1)

get_seurat_pcs <- function(counts_mat, npcs) {
  
  seurat_object <- CreateSeuratObject(counts = counts_mat)
  
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  
  seurat_object <- RunPCA(
    seurat_object, npcs = npcs, features = Features(seurat_object)
  )  
  
  # This is an orthogonal matrix but not orthonormal
  cell_embeddings <- seurat_object@reductions$pca@cell.embeddings
  
  cell_embeddings_normed <- sweep(
    cell_embeddings, 2, sqrt(colSums(cell_embeddings^2)),`/`
  )
  
  return(cell_embeddings_normed)
  
}

pbmc_cell_embeddings <- get_seurat_pcs(Matrix::t(counts), 10)

readr::write_rds(
  droplets_cell_embeddings,
  "~/Documents/data/fastglmpca/experiment_results/seurat_pca_pbmc.rds"
)

