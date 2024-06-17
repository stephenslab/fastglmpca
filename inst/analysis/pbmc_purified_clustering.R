# TO DO: Explain here what this script is for, and how to use it.
library(Matrix)
library(uwot)
library(ggplot2)
library(cowplot)
pbmc_colors <- c("dodgerblue","forestgreen","darkmagenta","salmon",
                 "gray","gold","yellow","orange","tomato","red")
set.seed(1)
load("results.RData")
k <- 10
V_scgbm      <- pbmc_purified_results$scGBM$fit$V
V_glmpca     <- pbmc_purified_results$glmpca$fit$V
V_fastglmpca <- pbmc_purified_results$fastglmpca_28_core$fit$V
colnames(V_scgbm) <- paste0("k",1:k)
colnames(V_glmpca) <- paste0("k",1:k)
colnames(V_fastglmpca) <- paste0("k",1:k)
cell_type <- sapply(strsplit(rownames(V_fastglmpca),"-"),"[[",3)
cell_type <- factor(cell_type)

# Plot the improvement in the solutions over time.
# TO DO.

# Plot the 10 PCs.
pc_plot <- function (V, k = 1:2, title = "") {
  k1 <- k[1]
  k2 <- k[2]
  pdat <- data.frame(cell_type = cell_type,
                     d1 = V[,k1],
                     d2 = V[,k2])
  return(ggplot(pdat,aes(x = d1,y = d2,color = cell_type)) +
         geom_point(size = 0.5,show.legend = FALSE) +
         scale_color_manual(values = pbmc_colors) +
         labs(x = paste("PC",k1),
              y = paste("PC",k2),
              title = title) +
         theme_cowplot(font_size = 9) +
         theme(plot.title = element_text(face = "plain",size = 9)))
}

p1 <- pc_plot(V_scgbm,1:2,title = "scGBM")
p2 <- pc_plot(V_scgbm,3:4,title = "scGBM")
p3 <- pc_plot(V_scgbm,5:6,title = "scGBM")
p4 <- pc_plot(V_scgbm,7:8,title = "scGBM")
p5 <- pc_plot(V_scgbm,9:10,title = "scGBM")

p6 <- pc_plot(V_glmpca,1:2,title = "glmpca")
p7 <- pc_plot(V_glmpca,3:4,title = "glmpca")
p8 <- pc_plot(V_glmpca,5:6,title = "glmpca")
p9 <- pc_plot(V_glmpca,7:8,title = "glmpca")
p10 <- pc_plot(V_glmpca,9:10,title = "glmpca")

p11 <- pc_plot(V_fastglmpca,1:2,title = "fastglmpca")
p12 <- pc_plot(V_fastglmpca,3:4,title = "fastglmpca")
p13 <- pc_plot(V_fastglmpca,5:6,title = "fastglmpca")
p14 <- pc_plot(V_fastglmpca,7:8,title = "fastglmpca")
p15 <- pc_plot(V_fastglmpca,9:10,title = "fastglmpca")

ggsave("pbmc_purified_pcs_k10.png",
       plot_grid(p1,p2,p3,p4,p5,
                 p6,p7,p8,p9,p10,
                 p11,p12,p13,p14,p15,
                 byrow = FALSE,nrow = 5,ncol = 3),
       height = 10,width = 7,dpi = 500,bg = "white")
       
# Project the fastglmpca PCs into 2-d using umap.
umap_plot <- function (Y, clusters, colors, title = "") {
  pdat <- data.frame(cluster = clusters,
                     umap1   = Y[,1],
                     umap2   = Y[,2])
  return(ggplot(pdat,aes(x = umap1,y = umap2,color = cluster)) +
         geom_point(size = 0.5) +
         scale_color_manual(values = colors) +
         ggtitle(title) +
         theme_cowplot(font_size = 9) +
         theme(plot.title = element_text(face = "plain",size = 9)))
         
}
res <- umap(V_fastglmpca,n_neighbors = 20,n_threads = 4,verbose = TRUE)
colnames(res) <- c("umap1","umap2")
p1 <- umap_plot(res,cell_type,pbmc_colors,title = "FACS")
