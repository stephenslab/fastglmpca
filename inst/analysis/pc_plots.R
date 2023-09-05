library(Matrix)
library(ggplot2)
library(cowplot)

load("results.RData")

# Compare pbmc, K = 2, glmpca vs. fastglmpca.
pbmc_colors <- c("forestgreen", # CD14+
                 "dodgerblue",  # B cells
                 "limegreen",   # CD34+
                 "gray",        # NK cells
                 "tomato",      # cytotoxic T cells
                 "darkmagenta", # dendritic
                 "gold")        # T cells
#load("~/git/fastTopics-experiments/data/pbmc_68k.RData")
load("~/Downloads/pbmc_68k.RData")
samples$celltype <- as.factor(dplyr::if_else(
  samples$celltype %in% c(
    "CD4+ T Helper2","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T",
    "CD4+/CD45RO+ Memory","CD8+/CD45RA+ Naive Cytotoxic"
  ),
  "T Cell",
  samples$celltype
))
rm(counts)
res1 <- pbmc_res_list$glmpca$`2_factors`
res2 <- pbmc_res_list$fastglmpca_fit_28_cores$`2_factors`
res3 <- pbmc_res_list$scGBM$`2_factors`
plot(res1$V,res2$V,pch = 20,cex = 0.75,
     xlab = "glmpca V",ylab = "fastglmpca V",
     col = "dodgerblue")
abline(a = 0,b = 1,col = "black",lty = "dashed")
pdat1 <- data.frame(celltype = samples$celltype,
                    PC1 = res1$V[,1],
                    PC2 = res1$V[,2])
pdat2 <- data.frame(celltype = samples$celltype,
                    PC1 = res2$V[,1],
                    PC2 = res2$V[,2])
pdat3 <- data.frame(celltype = samples$celltype,
                    PC1 = -res3$V[,2],
                    PC2 = -res3$V[,1])
p1 <- ggplot(pdat1,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("glmpca") + 
  theme_cowplot(font_size = 10)
p2 <- ggplot(pdat2,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("fastglmpca") + 
  theme_cowplot(font_size = 10)
p3 <- ggplot(pdat3,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("scGBM") + 
  theme_cowplot(font_size = 10)

plot_grid(p1,p2,p3,nrow = 2,ncol = 2)

# Compare droplet, K = 2, glmpca vs. fastglmpca.
tissue_colors <- c("royalblue",   # basal
                   "firebrick",   # ciliated
                   "forestgreen", # club
                   "gold",        # goblet
                   "darkmagenta", # ionocyte
                   "darkorange",  # neuroendocrine
                   "skyblue")     # tuft
#load("~/git/fastTopics-experiments/data/droplet.RData")
load("~/Downloads/droplet.RData")
res1 <- droplets_res_list$glmpca$`2_factors`
res2 <- droplets_res_list$fastglmpca_fit_28_cores$`2_factors`
res3 <- droplets_res_list$scGBM$`2_factors`
plot(res1$V,c(-res2$V[,2],res2$V[,1]),pch = 20,cex = 0.75,
     xlab = "glmpca V",ylab = "fastglmpca V",col = "dodgerblue")
abline(a = 0,b = 1,col = "black",lty = "dashed")
pdat1 <- data.frame(tissue = samples$tissue,
                    PC1 = res1$V[,1],
                    PC2 = res1$V[,2])
pdat2 <- data.frame(tissue = samples$tissue,
                    PC1 = -res2$V[,2],
                    PC2 = res2$V[,1])
pdat3 <- data.frame(tissue = samples$tissue,
                    PC1 = res3$V[,1],
                    PC2 = res3$V[,2])
p1 <- ggplot(pdat1,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("glmpca") + 
  theme_cowplot(font_size = 10)
p2 <- ggplot(pdat2,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("fastglmpca") + 
  theme_cowplot(font_size = 10)
p3 <- ggplot(pdat3,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("scGBM") +
  xlim(c(-0.05,0.03)) +
  ylim(c(-0.025,0.02)) +
  theme_cowplot(font_size = 10)
plot_grid(p1,p2,p3,nrow = 2,ncol = 2)

# Compare pbmc, K = 3, glmpca vs. fastglmpca.
# res1 <- readRDS(paste("pbmc_glmpca_fit_3_factors_avagrad_optimizer",
#                       "minibatch_stochastic_10_hrs.rds",sep="_"))
# res1 <- fastglmpca:::orthonormalize(as.matrix(res1$loadings),
#                                     as.matrix(res1$factors))
# res2 <- readRDS("pbmc_fastglmpca_fit_1_core_3_factors_10_hrs.rds")
# res2 <- fastglmpca:::orthonormalize(res2$U[,3:5],res2$V[,3:5])
# plot(res1$V,-res2$V[,c(3,1,2)],pch = 20)
# abline(a = 0,b = 1,col = "cyan",lty = "dotted")
