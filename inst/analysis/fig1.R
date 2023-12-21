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
  ggtitle("PBMC 68k glmpca") +
  theme_cowplot(font_size = 10)
p2 <- ggplot(pdat2,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k fastglmpca") +
  theme_cowplot(font_size = 10)
p3 <- ggplot(pdat3,aes(x = PC1,y = PC2,color = celltype)) +
  geom_point(size = 1) +
  scale_color_manual(values = pbmc_colors) +
  ggtitle("PBMC 68k scGBM") +
  theme_cowplot(font_size = 10)

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
p4 <- ggplot(pdat1,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca") +
  theme_cowplot(font_size = 10)
p5 <- ggplot(pdat2,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca") +
  theme_cowplot(font_size = 10)
p6 <- ggplot(pdat3,aes(x = PC1,y = PC2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets scGBM") +
  xlim(c(-0.05,0.03)) +
  ylim(c(-0.025,0.02)) +
  theme_cowplot(font_size = 10)
#plot_grid(p1,p2,p3,nrow = 2,ncol = 2)
library(ggpubr)
agg1 <- ggarrange(
  p1, p2, p3, nrow = 1, ncol = 3, common.legend = TRUE, legend = "right"
)
agg2 <- ggarrange(
  p4, p5, p6, nrow = 1, ncol = 3, common.legend = TRUE, legend = "right"
)


create_dataset_df <- function(
    results,
    factors
) {

  lik_df <- data.frame()

  for (n_factor in factors) {

    factors_str <- glue::glue("{n_factor}_factors")

    fastglmpca_28cores_fit <- results[["fastglmpca_fit_28_cores"]][[factors_str]]

    fastglmpca_1core_fit <- results[["fastglmpca_fit_1_core"]][[factors_str]]

    glmpca_fit <- results[["glmpca"]][[factors_str]]

    scGBM_fit <- results[["scGBM"]][[factors_str]]

    lik_factor_df <- data.frame(
      factors = rep(
        n_factor,
        length(fastglmpca_28cores_fit$time) + length(fastglmpca_1core_fit$time) +
          length(scGBM_fit$time) + length(glmpca_fit$time)
      ),
      Algorithm = c(
        rep("fastglmpca-28core", length(fastglmpca_28cores_fit$time)),
        rep("fastglmpca-1core", length(fastglmpca_1core_fit$time)),
        rep("scGBM", length(scGBM_fit$time)),
        rep("glmpca-avagrad-sgd", length(glmpca_fit$time))
      ),
      loglik = c(
        fastglmpca_28cores_fit$loglik,
        fastglmpca_1core_fit$loglik,
        scGBM_fit$loglik,
        glmpca_fit$loglik
      ),
      time = c(
        fastglmpca_28cores_fit$time,
        fastglmpca_1core_fit$time,
        scGBM_fit$time,
        glmpca_fit$time
      )
    )

    best_ll <- max(lik_factor_df$loglik)

    lik_factor_df <- lik_factor_df %>%
      dplyr::mutate(dist_from_best_ll = abs(loglik - best_ll) + .01)

    lik_df <- rbind(lik_df, lik_factor_df)

  }

  return(lik_df)

}

create_plot_list <- function(
    dataset_df,
    dataset,
    factors
) {

  l10_list <- list()
  #std_list <- list()

  for (n_factor in factors) {

    factor_df <- dataset_df %>%
      dplyr::filter(factors == n_factor)

    factor_df$Algorithm <- factor(
      factor_df$Algorithm,
      levels=c(
        "scGBM", "glmpca-avagrad-sgd", "fastglmpca-1core", "fastglmpca-28core"
        )
    )

    g_l10 <- ggplot(data = factor_df) +
      geom_line(aes(x = time / 3600, y = dist_from_best_ll, color = Algorithm)) +
      scale_y_continuous(trans = "log10") +
      xlab("Time (hours)") +
      ylab("Distance from Best Log-Likelihood") +
      ggtitle(glue::glue("{dataset} {n_factor} Factors")) +
      scale_color_manual(
        values = c(
          "royalblue",
          "darkorange",
          "forestgreen",
          "gold"
        )
      ) +
      cowplot::theme_cowplot()

    l10_list[[glue::glue("l10_{n_factor}")]] <- g_l10

  }

  return(
    c(
      l10_list
      #std_list
    )
  )

}

droplet_df <- create_dataset_df(
  droplets_res_list,
  c(2)
)

droplet_plots_list <- create_plot_list(
  droplet_df,
  "Droplets",
  c(2)
)

pbmc_df <- create_dataset_df(
  pbmc_res_list,
  c(2)
)

pbmc_plots_list <- create_plot_list(
  pbmc_df,
  "PBMC 68k",
  c(2)
)

agg3 <- ggarrange(
  droplet_plots_list$l10_2, pbmc_plots_list$l10_2, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right"
)

ggarrange(agg1, agg2, agg3, nrow = 3, ncol = 1)

# Compare pbmc, K = 3, glmpca vs. fastglmpca.
# res1 <- readRDS(paste("pbmc_glmpca_fit_3_factors_avagrad_optimizer",
#                       "minibatch_stochastic_10_hrs.rds",sep="_"))
# res1 <- fastglmpca:::orthonormalize(as.matrix(res1$loadings),
#                                     as.matrix(res1$factors))
# res2 <- readRDS("pbmc_fastglmpca_fit_1_core_3_factors_10_hrs.rds")
# res2 <- fastglmpca:::orthonormalize(res2$U[,3:5],res2$V[,3:5])
# plot(res1$V,-res2$V[,c(3,1,2)],pch = 20)
# abline(a = 0,b = 1,col = "cyan",lty = "dotted")
