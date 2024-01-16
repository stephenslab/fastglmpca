# Here, I want to see if I can find anything interesting in the
# other datasets for different values of K. I'm curious if 
# fastglmpca is much better for K = 4, for instance
# this could be for either droplets of pbmc

fastglmpca_10_factors <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/droplets_fastglmpca_fit_10_factors_5105_iter_28_cores_dec_23.rds"
)

glmpca_10_factors <- readr::read_rds(
  "~/Documents/data/fastglmpca/experiment_results/droplets_glmpca_fit_10_factors_10_hrs_avagrad_optimizer_minibatch_stochastic_dec_23.rds"
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
  ggtitle("Droplets fastglmpca PCs 1 and 2") +
  cowplot::theme_cowplot(font_size = 10)
g2 <- ggplot(glmpca_df,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 1 and 2") +
  cowplot::theme_cowplot(font_size = 10)

g3 <- ggplot(fastglmpca_df,aes(x = PC3,y = PC4,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 3 and 4") +
  cowplot::theme_cowplot(font_size = 10)
g4 <- ggplot(glmpca_df,aes(x = PC3,y = PC4,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 3 and 4") +
  cowplot::theme_cowplot(font_size = 10)

g5 <- ggplot(fastglmpca_df,aes(x = PC5,y = PC6,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 5 and 6") +
  cowplot::theme_cowplot(font_size = 10)
g6 <- ggplot(glmpca_df,aes(x = PC5,y = PC6,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 5 and 6") +
  cowplot::theme_cowplot(font_size = 10)

g7 <- ggplot(fastglmpca_df,aes(x = PC7,y = PC8,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 7 and 8") +
  cowplot::theme_cowplot(font_size = 10)
g8 <- ggplot(glmpca_df,aes(x = PC7,y = PC8,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 7 and 8") +
  cowplot::theme_cowplot(font_size = 10)

g9 <- ggplot(fastglmpca_df,aes(x = PC9,y = PC10,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 9 and 10") +
  cowplot::theme_cowplot(font_size = 10)
g10 <- ggplot(glmpca_df,aes(x = PC9,y = PC10,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 9 and 10") +
  cowplot::theme_cowplot(font_size = 10)

library(ggpubr)
ggarrange(
  g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, nrow = 5, ncol = 2, common.legend = TRUE, legend = "right"
  )




