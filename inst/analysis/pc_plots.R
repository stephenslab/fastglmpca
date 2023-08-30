library(Matrix)
library(ggplot2)
library(cowplot)

# Compare pbmc, K = 2, glmpca vs. fastglmpca.
cell_type_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
                      "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
                      "#ffff99")
#load("~/git/fastTopics-experiments/data/pbmc_68k.RData")
load("~/Downloads/pbmc_68k.RData")
rm(counts)
res1 <- readRDS(paste("pbmc_glmpca_fit_2_factors_avagrad_optimizer",
                      "minibatch_stochastic_10_hrs.rds",sep="_"))
res1 <- fastglmpca:::orthonormalize(as.matrix(res1$loadings),
                                    as.matrix(res1$factors))
res2 <- readRDS("pbmc_fastglmpca_fit_28_cores_2_factors_10_hrs.rds")
res2 <- fastglmpca:::orthonormalize(res2$U[,3:4],res2$V[,3:4])
res3 <- readRDS("pbmc_scGBM_fit_2_factors_no_beta_infer_10_hrs.rds")
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
  scale_color_manual(values = cell_type_colors) +
  ggtitle("glmpca") + 
  theme_cowplot(font_size = 10)
p2 <- ggplot(pdat2,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_type_colors) +
  ggtitle("fastglmpca") + 
  theme_cowplot(font_size = 10)
p3 <- ggplot(pdat3,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = cell_type_colors) +
  ggtitle("scGBM") + 
  theme_cowplot(font_size = 10)

plot_grid(p1,p2,p3,nrow = 2,ncol = 2)

# Compare droplet, K = 2, glmpca vs. fastglmpca.
tissue_colors <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00",
                   "#ffff33","#a65628")
#load("~/git/fastTopics-experiments/data/droplet.RData")
load("~/Downloads/droplet.RData")
res1 <- readRDS(paste("droplets_glmpca_fit_2_factors_avagrad_optimizer",
                      "minibatch_stochastic_10_hrs.rds",sep="_"))
res1 <- fastglmpca:::orthonormalize(as.matrix(res1$loadings),
                                    as.matrix(res1$factors))
res2 <- readRDS("droplets_fastglmpca_fit_28_cores_2_factors_10_hrs.rds")
res2 <- fastglmpca:::orthonormalize(res2$U[,3:4],res2$V[,3:4])
res3 <- readRDS("droplets_scGBM_fit_2_factors_no_beta_infer_10_hrs.rds")
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
