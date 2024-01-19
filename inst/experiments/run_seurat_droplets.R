load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

library(Seurat)
set.seed(1)

get_seurat_2pcs <- function(counts_mat) {
  
  seurat_object <- CreateSeuratObject(counts = counts_mat)
  
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- ScaleData(seurat_object)
    
  seurat_object <- RunPCA(
    seurat_object, npcs = 2, features = Features(seurat_object)
  )  
  
  # This is an orthogonal matrix but not orthonormal
  cell_embeddings <- seurat_object@reductions$pca@cell.embeddings
  
  cell_embeddings_normed <- sweep(
    cell_embeddings, 2, sqrt(colSums(cell_embeddings^2)),`/`
  )
  
  return(cell_embeddings_normed)
  
}

droplets_cell_embeddings <- get_seurat_2pcs(Matrix::t(counts))

readr::write_rds(
  droplets_cell_embeddings,
  "~/Documents/data/fastglmpca/experiment_results/seurat_pca_droplets.rds"
)

tissue_colors <- c("royalblue",   # basal
                   "firebrick",   # ciliated
                   "forestgreen", # club
                   "gold",        # goblet
                   "darkmagenta", # ionocyte
                   "darkorange",  # neuroendocrine
                   "skyblue")     # tuft

pdat1 <- data.frame(tissue = samples$tissue,
                    PC1 = droplets_cell_embeddings[,1],
                    PC2 = droplets_cell_embeddings[,2])

# Create a ggplot2 scatter plot
library(ggplot2)
droplets_pca_plot <- ggplot(pdat1,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets Seurat") +
  cowplot::theme_cowplot(font_size = 10)
