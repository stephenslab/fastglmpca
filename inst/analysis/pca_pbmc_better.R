library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
pbmc_colors <- c("forestgreen", # CD14+
                 "dodgerblue",  # B cells
                 "limegreen",   # CD34+
                 "gray",        # NK cells
                 "tomato",      # cytotoxic T cells
                 "darkmagenta", # dendritic
                 "gold")        # T cells
load("~/git/fastTopics-experiments/data/pbmc_68k.RData")
rm(counts)
load("results.RData")
V <- pbmc_res_list$fastglmpca_28_cores$PCs_10$V
pdat <- cbind(samples["celltype"],V[,1:2])
names(pdat) <- c("celltype","PC1","PC2")
pdat$celltype <- as.factor(if_else(
  pdat$celltype %in%
    c("CD4+ T Helper2","CD4+/CD25 T Reg","CD4+/CD45RA+/CD25- Naive T",
      "CD4+/CD45RO+ Memory","CD8+/CD45RA+ Naive Cytotoxic"),
    "T Cell",
    pdat$celltype))
p1 <- ggplot(pdat,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  theme_cowplot(font_size = 10)
p2 <- ggplot(pdat,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  xlim(-0.004,0.013) +
  ylim(-0.0115,0.01) +
  scale_color_manual(values = pbmc_colors) +
  theme_cowplot(font_size = 10)
ggsave("pbmc_fastglmpca_k10.pdf",plot_grid(p1,p2),height = 2.5,width = 8)
