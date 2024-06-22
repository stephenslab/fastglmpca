# Generate plots summarizing results on the "purified" PBMC data, with
# K = 10.
library(Matrix)
library(uwot)
library(ggplot2)
library(cowplot)
methods_colors <- c("darkorange","magenta","dodgerblue","darkblue")
pbmc_colors    <- c("dodgerblue","forestgreen","darkmagenta","gray",
                    "gold","tomato")
cluster_colors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
                    "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a")
set.seed(1)
load("results.RData")

# Plot the improvement in the log-likelihoods over time.
n_glmpca <- length(pbmc_purified_results$glmpca$fit$loglik)
pdat <-
  rbind(
      data.frame(method = "scGBM",
                   time   = pbmc_purified_results$scGBM$fit$time,
                   loglik = pbmc_purified_results$scGBM$fit$loglik),
        data.frame(method = "glmpca",
                   time   = seq(0,10*60^2,length.out = n_glmpca),
                   loglik = pbmc_purified_results$glmpca$fit$loglik),
        data.frame(method = "fastglmpca_1core",
                   time   = cumsum(pbmc_purified_results$fastglmpca_1_core$fit$progress$time),
                   loglik = pbmc_purified_results$fastglmpca_1_core$fit$progress$loglik),
        data.frame(method = "fastglmpca_28core",
                   time   = cumsum(pbmc_purified_results$fastglmpca_28_core$fit$progress$time),
                   loglik = pbmc_purified_results$fastglmpca_28_core$fit$progress$loglik))
pdat <- transform(pdat,
                  method  = factor(method),
                  time    = time/60^2,
                  loglik  = max(loglik) - loglik + 1)
p1 <- ggplot(pdat,aes(x = time,y = loglik,color = method)) +
  geom_line(size = 0.7) +
  scale_x_continuous(breaks = seq(0,10)) +
  scale_y_continuous(trans = "log10",breaks = 10^seq(0,8)) +
  scale_color_manual(values = methods_colors) +
  labs(x = "running time (h)",y = "distance from best loglik") +
  theme_cowplot(font_size = 10)
ggsave("pbmc_purified_loglik.pdf",p1,height = 3,width = 4.6)

# Plot the change in the clustering (NMI, ARI) over time.
n_glmpca <- length(pbmc_purified_results$glmpca$cluster_metrics_by_iter$nmi)
n1 <- length(pbmc_purified_results$fastglmpca_1_core$fit$progress$time)
n2 <- length(pbmc_purified_results$fastglmpca_28_core$fit$progress$time)
pdat <-
  rbind(
    data.frame(method = "scGBM",
               time = c(0,pbmc_purified_results$scGBM$fit$time),
               nmi  = pbmc_purified_results$scGBM$cluster_metrics_by_iter$nmi,
               ari  = pbmc_purified_results$scGBM$cluster_metrics_by_iter$ari),
    data.frame(method = "glmpca",
               time = seq(0,10*60^2,length.out = n_glmpca),
               nmi  = pbmc_purified_results$glmpca$cluster_metrics_by_iter$nmi,
               ari  = pbmc_purified_results$glmpca$cluster_metrics_by_iter$ari),
    data.frame(method = "fastglmpca_1core",
               time = cumsum(pbmc_purified_results$fastglmpca_1_core$fit$progress$time)[-1],
               nmi  = pbmc_purified_results$fastglmpca_1_core$cluster_metrics_by_iter$nmi,
               ari  = pbmc_purified_results$fastglmpca_1_core$cluster_metrics_by_iter$ari),
    data.frame(method = "fastglmpca_28core",
               time = cumsum(pbmc_purified_results$fastglmpca_28_core$fit$progress$time)[-1],
               nmi  = pbmc_purified_results$fastglmpca_28_core$cluster_metrics_by_iter$nmi,
               ari  = pbmc_purified_results$fastglmpca_28_core$cluster_metrics_by_iter$ari))
pdat <- transform(pdat,
                  method = factor(method),
                  time   = time/60^2)
p1 <- ggplot(pdat,aes(x = time,y = nmi,color = method)) +
  geom_line(size = 0.7) +
  scale_x_continuous(breaks = seq(0,10)) +
  scale_color_manual(values = methods_colors) +
  labs(x = "running time (h)",y = "NMI") +
  theme_cowplot(font_size = 10)
p2 <- ggplot(pdat,aes(x = time,y = ari,color = method)) +
  geom_line(size = 0.7) +
  scale_x_continuous(breaks = seq(0,10)) +
  scale_color_manual(values = methods_colors) +
  labs(x = "running time (h)",y = "ARI") +
  theme_cowplot(font_size = 10)
ggsave("pbmc_purified_nmi_ari.pdf",
       plot_grid(p1,p2),
       height = 2.1,width = 7)

# Plot the 10 PCs.
pc_plot <- function (V, k = 1:2, title = "", show.legend = FALSE) {
  k1 <- k[1]
  k2 <- k[2]
  pdat <- data.frame(cell_type = cell_type,
                     d1 = V[,k1],
                     d2 = V[,k2])
  return(ggplot(pdat,aes(x = d1,y = d2,color = cell_type)) +
         geom_point(size = 0.5,show.legend = show.legend) +
         scale_color_manual(values = pbmc_colors) +
         labs(x = paste("PC",k1),
              y = paste("PC",k2),
              title = title) +
         theme_cowplot(font_size = 9) +
         theme(plot.title = element_text(face = "plain",size = 9)))
}

k <- 10
V_scgbm      <- pbmc_purified_results$scGBM$fit$V
V_glmpca     <- pbmc_purified_results$glmpca$fit$V
V_fastglmpca <- pbmc_purified_results$fastglmpca_28_core$fit$V
colnames(V_scgbm)      <- paste0("k",1:k)
colnames(V_glmpca)     <- paste0("k",1:k)
colnames(V_fastglmpca) <- paste0("k",1:k)
cell_type <- sapply(strsplit(rownames(V_fastglmpca),"-"),"[[",3)
cell_type[cell_type == "memory_t" |
          cell_type == "regulatory_t" |
          cell_type == "cd4_t_helper" |
          cell_type == "naive_cytotoxic" |
          cell_type == "naive_t"] <- "t_cells"
cell_type <- factor(cell_type)

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

p16 <- pc_plot(V_fastglmpca,1:2,show.legend = TRUE)

ggsave("pbmc_purified_pcs_k10.png",
       plot_grid(p1,p2,p3,p4,p5,
                 p6,p7,p8,p9,p10,
                 p11,p12,p13,p14,p15,
                 byrow = FALSE,nrow = 5,ncol = 3),
       height = 10,width = 7,dpi = 500,bg = "white")
ggsave("pbmc_purified_pcs_k10_legend.png",p16,
       height = 3,width = 3,dpi = 500,bg = "white")
