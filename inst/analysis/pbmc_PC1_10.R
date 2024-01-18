fastglmpca_10_factors <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/pbmc_fastglmpca_fit_10_factors_358_iter_28_cores_dec_23.rds"
)

glmpca_10_factors <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/pbmc_glmpca_fit_10_factors_10_hrs_avagrad_optimizer_minibatch_stochastic_dec_23.rds"
)

scGMB_10_factors <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/pbmc_scGBM_fit_10_factors_no_beta_infer_10_hrs.rds"
)

glmpca_mod <- fastglmpca:::orthonormalize(
  as.matrix(glmpca_10_factors$loadings), as.matrix(glmpca_10_factors$factors)
)

load("~/Documents/data/fastglmpca/raw_data/pbmc_68k.RData")

samples$celltype <- as.factor(dplyr::if_else(
  samples$celltype %in% c(
    "CD4+ T Helper2","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T",
    "CD4+/CD45RO+ Memory","CD8+/CD45RA+ Naive Cytotoxic"
  ),
  "T Cell",
  samples$celltype
))

fastglmpca_df <- data.frame(
  celltype = samples$celltype,
  PC1 = fastglmpca_10_factors$V[,1],
  PC2 = fastglmpca_10_factors$V[,2],
  PC3 = fastglmpca_10_factors$V[,3],
  PC4 = fastglmpca_10_factors$V[,4],
  PC5 = fastglmpca_10_factors$V[,5],
  PC6 = fastglmpca_10_factors$V[,6],
  PC7 = fastglmpca_10_factors$V[,7],
  PC8 = fastglmpca_10_factors$V[,8],
  PC9 = fastglmpca_10_factors$V[,9],
  PC10 = fastglmpca_10_factors$V[,10]
)

glmpca_df <- data.frame(
  celltype = samples$celltype,
  PC1 = glmpca_mod$V[,1],
  PC2 = glmpca_mod$V[,2],
  PC3 = glmpca_mod$V[,3],
  PC4 = glmpca_mod$V[,4],
  PC5 = glmpca_mod$V[,5],
  PC6 = glmpca_mod$V[,6],
  PC7 = glmpca_mod$V[,7],
  PC8 = glmpca_mod$V[,8],
  PC9 = glmpca_mod$V[,9],
  PC10 = glmpca_mod$V[,10]
)

scGBM_df <- data.frame(
  celltype = samples$celltype,
  PC1 = scGMB_10_factors$V[,1],
  PC2 = scGMB_10_factors$V[,2],
  PC3 = scGMB_10_factors$V[,3],
  PC4 = scGMB_10_factors$V[,4],
  PC5 = scGMB_10_factors$V[,5],
  PC6 = scGMB_10_factors$V[,6],
  PC7 = scGMB_10_factors$V[,7],
  PC8 = scGMB_10_factors$V[,8],
  PC9 = scGMB_10_factors$V[,9],
  PC10 = scGMB_10_factors$V[,10]
)

library(ggplot2)

pbmc_colors <- c("forestgreen", # CD14+
                 "dodgerblue",  # B cells
                 "limegreen",   # CD34+
                 "gray",        # NK cells
                 "tomato",      # cytotoxic T cells
                 "darkmagenta", # dendritic
                 "gold")        # T cells


g1 <- ggplot(fastglmpca_df,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k fastglmpca PCs 1 & 2") +
  cowplot::theme_cowplot(font_size = 10)
g2 <- ggplot(glmpca_df,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k glmpca PCs 1 & 2") +
  cowplot::theme_cowplot(font_size = 10)
g3 <- ggplot(scGBM_df,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k scGBM PCs 1 & 2") +
  cowplot::theme_cowplot(font_size = 10)

g4 <- ggplot(fastglmpca_df,aes(x = PC3,y = PC4,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k fastglmpca PCs 3 & 4") +
  cowplot::theme_cowplot(font_size = 10)
g5 <- ggplot(glmpca_df,aes(x = PC3,y = PC4,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k glmpca PCs 3 & 4") +
  cowplot::theme_cowplot(font_size = 10)
g6 <- ggplot(scGBM_df,aes(x = PC3,y = PC4,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k scGBM PCs 3 & 4") +
  cowplot::theme_cowplot(font_size = 10)

g7 <- ggplot(fastglmpca_df,aes(x = PC5,y = PC6,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k fastglmpca PCs 5 & 6") +
  cowplot::theme_cowplot(font_size = 10)
g8 <- ggplot(glmpca_df,aes(x = PC5,y = PC6,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k glmpca PCs 5 & 6") +
  cowplot::theme_cowplot(font_size = 10)
g9 <- ggplot(scGBM_df,aes(x = PC5,y = PC6,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k scGBM PCs 5 & 6") +
  cowplot::theme_cowplot(font_size = 10)

g10 <- ggplot(fastglmpca_df,aes(x = PC7,y = PC8,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k fastglmpca PCs 7 & 8") +
  cowplot::theme_cowplot(font_size = 10)
g11 <- ggplot(glmpca_df,aes(x = PC7,y = PC8,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k glmpca PCs 7 & 8") +
  cowplot::theme_cowplot(font_size = 10)
g12 <- ggplot(scGBM_df,aes(x = PC7,y = PC8,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k scGBM PCs 7 & 8") +
  cowplot::theme_cowplot(font_size = 10)

g13 <- ggplot(fastglmpca_df,aes(x = PC9,y = PC10,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k fastglmpca PCs 9 & 10") +
  cowplot::theme_cowplot(font_size = 10)
g14 <- ggplot(glmpca_df,aes(x = PC9,y = PC10,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k glmpca PCs 9 & 10") +
  cowplot::theme_cowplot(font_size = 10)
g15 <- ggplot(scGBM_df,aes(x = PC9,y = PC10,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k scGBM PCs 9 & 10") +
  cowplot::theme_cowplot(font_size = 10)

png(
  file="~/Documents/fastglmpca/inst/scratch/pbmc_PC1_10.png", 
  height = 100 * 11,
  width = 100 * 8.5
  )

library(ggpubr)
ggarrange(
  g1, g2, g3, 
  g4, g5, g6, 
  g7, g8, g9, 
  g10, g11, g12,
  g13, g14, g15,
  nrow = 5, ncol = 3, 
  common.legend = TRUE, legend = "right"
)
dev.off()

