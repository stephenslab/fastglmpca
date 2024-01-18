library(Matrix)
library(ggplot2)
library(cowplot)
library(dplyr)

load("~/Documents/fastglmpca/inst/analysis/results.RData")

pbmc_cell_embeddings <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/seurat_pca_pbmc.rds"
)

pbmc_df <- data.frame(pbmc_cell_embeddings)

droplets_cell_embeddings <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/seurat_pca_droplets.rds"
)

droplets_df <- data.frame(droplets_cell_embeddings)

droplets_df <- droplets_df %>%
  dplyr::rename(PC1 = PC_1, PC2 = PC_2)

# Compare pbmc, K = 2, glmpca vs. fastglmpca.
pbmc_colors <- c("forestgreen", # CD14+
                 "dodgerblue",  # B cells
                 "limegreen",   # CD34+
                 "gray",        # NK cells
                 "tomato",      # cytotoxic T cells
                 "darkmagenta", # dendritic
                 "gold")        # T cells
#load("~/git/fastTopics-experiments/data/pbmc_68k.RData")
load("~/Documents/data/fastglmpca/raw_data/pbmc_68k.RData")
samples$celltype <- as.factor(dplyr::if_else(
  samples$celltype %in% c(
    "CD4+ T Helper2","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T",
    "CD4+/CD45RO+ Memory","CD8+/CD45RA+ Naive Cytotoxic"
  ),
  "T Cell",
  samples$celltype
))

pbmc_df$celltype <- samples$celltype

res1 <- pbmc_res_list$glmpca$`2_factors`
res2 <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/pbmc_fastglmpca_fit_2_factors_1921_iter_28_cores_dec_23.rds"
)
res3 <- pbmc_res_list$scGBM$`2_factors`

pdat1 <- data.frame(celltype = samples$celltype,
                    PC1 = res1$V[,1],
                    PC2 = res1$V[,2])
pdat2 <- data.frame(celltype = samples$celltype,
                    PC1 = res2$V[,1],
                    PC2 = res2$V[,2])
pdat3 <- data.frame(celltype = samples$celltype,
                    PC1 = -res3$V[,2],
                    PC2 = -res3$V[,1])

p1 <- ggplot(pbmc_df,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k Seurat") +
  cowplot::theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
p2 <- ggplot(pdat1,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k glmpca") +
  theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
p3 <- ggplot(pdat2,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k fastglmpca") +
  theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
p4 <- ggplot(pdat3,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k scGBM") +
  theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)

#plot_grid(p1,p2,p3,nrow = 2,ncol = 2)

# Compare droplet, K = 2, glmpca vs. fastglmpca.
tissue_colors <- c("royalblue",   # basal
                   "firebrick",   # ciliated
                   "forestgreen", # club
                   "gold",        # goblet
                   "darkmagenta", # ionocyte
                   "darkorange",  # neuroendocrine
                   "skyblue")     # tuft
#load("~/git/fastTopics-experiments/data/droplet.RData")
load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

res1 <- droplets_res_list$glmpca$`2_factors`
res2 <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/droplets_fastglmpca_fit_2_factors_23135_iter_28_cores_dec_23.rds"
)
res3 <- droplets_res_list$scGBM$`2_factors`

droplets_df$tissue <- samples$tissue

pdat1 <- data.frame(tissue = samples$tissue,
                    PC1 = res1$V[,1],
                    PC2 = res1$V[,2])
pdat2 <- data.frame(tissue = samples$tissue,
                    PC1 = -res2$V[,2],
                    PC2 = res2$V[,1])
pdat3 <- data.frame(tissue = samples$tissue,
                    PC1 = res3$V[,1],
                    PC2 = res3$V[,2])

p5 <- ggplot(droplets_df,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets Seurat") +
  cowplot::theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
p6 <- ggplot(pdat1,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca") +
  theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
p7 <- ggplot(pdat2,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca") +
  theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
p8 <- ggplot(pdat3,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets scGBM") +
  xlim(c(-0.05,0.03)) +
  ylim(c(-0.025,0.02)) +
  theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
#plot_grid(p1,p2,p3,nrow = 2,ncol = 2)
library(ggpubr)
agg1 <- ggarrange(
  p1, p2, p3, p4, nrow = 1, ncol = 4, 
  common.legend = TRUE, legend = "right",
  labels = c("A", "", "")
)
agg2 <- ggarrange(
  p5, p6, p7, p8, nrow = 1, ncol = 4, 
  common.legend = TRUE, legend = "right",
  labels = c("B", "", "")
)

ggarrange(agg1, agg2, nrow = 2, ncol = 1)
ggsave(
  "~/Documents/fastglmpca/inst/scratch/PCs_1_2.pdf",
  device = "pdf",
  width = 16,
  height = 8
  )

