# Now, want to make this for the pbmc dataset
load("~/Documents/data/fastglmpca/pbmc_68k.RData")

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

pbmc_cell_embeddings <- get_seurat_2pcs(
  Matrix::t(counts)
)

# pbmc_colors <- c("forestgreen", # CD14+
#                  "dodgerblue",  # B cells
#                  "limegreen",   # CD34+
#                  "gray",        # NK cells
#                  "tomato",      # cytotoxic T cells
#                  "darkmagenta", # dendritic
#                  "gold")        # T cells
# 
# samples$celltype <- as.factor(dplyr::if_else(
#   samples$celltype %in% c(
#     "CD4+ T Helper2","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T",
#     "CD4+/CD45RO+ Memory","CD8+/CD45RA+ Naive Cytotoxic"
#   ),
#   "T Cell",
#   samples$celltype
# ))

pdat <- data.frame(tissue = samples$celltype,
                    PC1 = pbmc_cell_embeddings[,1],
                    PC2 = pbmc_cell_embeddings[,2])

# pbmc_pca_plot <- ggplot(pdat2,aes(x = PC1,y = PC2,color = celltype)) +
#   geom_point(size = 1) +
#   scale_color_manual(values = pbmc_colors) +
#   ggtitle("PBMC 68k Seurat") +
#   cowplot::theme_cowplot(font_size = 10)

readr::write_rds(pdat, "seurat_pca_droplets.rds")

