library(Matrix)
library(uwot)
library(ggplot2)
library(cowplot)
pbmc_colors <- c("dodgerblue","forestgreen","darkmagenta","salmon",
                 "gray","gold","yellow","orange","tomato","red")
set.seed(1)
load("results.RData")
V <- pbmc_purified_results$fastglmpca_28_cores[["25_factors"]]$V
cell_type <- sapply(strsplit(rownames(V),"-"),"[[",3)
cell_type <- factor(cell_type)

# Plot the first 2 PCs.
pdat <- cbind(data.frame(cell_type = cell_type),V)
p1 <- ggplot(pdat,aes(x = k_1,y = k_2,color = cell_type)) +
  geom_point() +
  scale_color_manual(values = pbmc_colors) +
  theme_cowplot(font_size = 10)

n <- nrow(V)
i <- sample(n,2e4)
cell_type <- cell_type[i]

# Project the PCs into 2-d using umap.
res <- umap(V[i,],n_neighbors = 20,n_threads = 4,verbose = TRUE)
colnames(res) <- c("umap1","umap2")
pdat <- cbind(data.frame(cell_type = cell_type),res)
p2 <- ggplot(pdat,aes(x = umap1,y = umap2,color = cell_type)) +
  geom_point() +
  scale_color_manual(values = pbmc_colors) +
  theme_cowplot(font_size = 10)
