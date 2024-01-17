fastglmpca_10_factors <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/droplets_fastglmpca_fit_10_factors_4500_iter_28_cores_dec_23.rds"
)

glmpca_10_factors <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/droplets_glmpca_fit_10_factors_10_hrs_avagrad_optimizer_minibatch_stochastic_dec_23.rds"
)

scGMB_10_factors <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/droplets_scGBM_fit_10_factors_no_beta_infer_10_hrs.rds"
)

load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

fastglmpca_df <- data.frame(
  tissue = samples$tissue,
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

glmpca_mod <- fastglmpca:::orthonormalize(
  as.matrix(glmpca_10_factors$loadings), as.matrix(glmpca_10_factors$factors)
)

glmpca_df <- data.frame(
  tissue = samples$tissue,
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
  tissue = samples$tissue,
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

tissue_colors <- c("royalblue",   # basal
                   "firebrick",   # ciliated
                   "forestgreen", # club
                   "gold",        # goblet
                   "darkmagenta", # ionocyte
                   "darkorange",  # neuroendocrine
                   "skyblue")     # tuft


g1 <- ggplot(fastglmpca_df,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 1 & 2") +
  cowplot::theme_cowplot(font_size = 10)
g2 <- ggplot(glmpca_df,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 1 & 2") +
  cowplot::theme_cowplot(font_size = 10)
g3 <- ggplot(scGBM_df,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets scGBM PCs 1 & 2") +
  cowplot::theme_cowplot(font_size = 10)

g4 <- ggplot(fastglmpca_df,aes(x = PC3,y = PC4,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 3 & 4") +
  cowplot::theme_cowplot(font_size = 10)
g5 <- ggplot(glmpca_df,aes(x = PC3,y = PC4,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 3 & 4") +
  cowplot::theme_cowplot(font_size = 10)
g6 <- ggplot(scGBM_df,aes(x = PC3,y = PC4,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets scGBM PCs 3 & 4") +
  cowplot::theme_cowplot(font_size = 10)

g7 <- ggplot(fastglmpca_df,aes(x = PC5,y = PC6,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 5 & 6") +
  cowplot::theme_cowplot(font_size = 10)
g8 <- ggplot(glmpca_df,aes(x = PC5,y = PC6,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 5 & 6") +
  cowplot::theme_cowplot(font_size = 10)
g9 <- ggplot(scGBM_df,aes(x = PC5,y = PC6,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets scGBM PCs 5 & 6") +
  cowplot::theme_cowplot(font_size = 10)

g10 <- ggplot(fastglmpca_df,aes(x = PC7,y = PC8,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 7 & 8") +
  cowplot::theme_cowplot(font_size = 10)
g11 <- ggplot(glmpca_df,aes(x = PC7,y = PC8,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 7 & 8") +
  cowplot::theme_cowplot(font_size = 10)
g12 <- ggplot(scGBM_df,aes(x = PC7,y = PC8,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets scGBM PCs 7 & 8") +
  cowplot::theme_cowplot(font_size = 10)

g13 <- ggplot(fastglmpca_df,aes(x = PC9,y = PC10,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 9 & 10") +
  cowplot::theme_cowplot(font_size = 10)
g14 <- ggplot(glmpca_df,aes(x = PC9,y = PC10,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 9 & 10") +
  cowplot::theme_cowplot(font_size = 10)
g15 <- ggplot(scGBM_df,aes(x = PC9,y = PC10,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets scGBM PCs 9 & 10") +
  cowplot::theme_cowplot(font_size = 10)

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

