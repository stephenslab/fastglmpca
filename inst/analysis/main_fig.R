load("~/Documents/fastglmpca/inst/analysis/results.RData")

fastglmpca_10_factors <- droplets_res_list$fastglmpca_28_cores$PCs_10

glmpca_mod <- droplets_res_list$glmpca$`10_factors`

scGMB_10_factors <- droplets_res_list$scGBM$`10_factors`

load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

fastglmpca_df <- data.frame(
  tissue = samples$tissue,
  PC3 = fastglmpca_10_factors$V[,3],
  PC4 = fastglmpca_10_factors$V[,4],
  PC5 = fastglmpca_10_factors$V[,5],
  PC6 = fastglmpca_10_factors$V[,6]
)

glmpca_df <- data.frame(
  tissue = samples$tissue,
  PC3 = glmpca_mod$V[,3],
  PC4 = glmpca_mod$V[,4],
  PC5 = glmpca_mod$V[,5],
  PC6 = glmpca_mod$V[,6]
)

scGBM_df <- data.frame(
  tissue = samples$tissue,
  PC3 = scGMB_10_factors$V[,3],
  PC4 = scGMB_10_factors$V[,4],
  PC5 = scGMB_10_factors$V[,5],
  PC6 = scGMB_10_factors$V[,6]
)

library(ggplot2)

tissue_colors <- c("royalblue",   # basal
                   "firebrick",   # ciliated
                   "forestgreen", # club
                   "gold",        # goblet
                   "darkmagenta", # ionocyte
                   "darkorange",  # neuroendocrine
                   "skyblue")     # tuft

g1 <- ggplot(fastglmpca_df,aes(x = PC3,y = PC4,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 3 & 4") +
  cowplot::theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
g2 <- ggplot(glmpca_df,aes(x = PC3,y = PC4,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 3 & 4") +
  cowplot::theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
g3 <- ggplot(scGBM_df,aes(x = PC3,y = PC4,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets scGBM PCs 3 & 4") +
  cowplot::theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)

g4 <- ggplot(fastglmpca_df,aes(x = PC5,y = PC6,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets fastglmpca PCs 5 & 6") +
  cowplot::theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
g5 <- ggplot(glmpca_df,aes(x = PC5,y = PC6,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets glmpca PCs 5 & 6") +
  cowplot::theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)
g6 <- ggplot(scGBM_df,aes(x = PC5,y = PC6,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("Droplets scGBM PCs 5 & 6") +
  cowplot::theme_cowplot(font_size = 10) +
  theme(aspect.ratio=1)

library(ggpubr)
pcs <- ggarrange(
  g1, g2, g3, 
  g4, g5, g6, 
  nrow = 2, ncol = 3, 
  common.legend = TRUE, legend = "right",
  labels = c("A", "", "", "B", "", "")
)

create_dataset_df <- function(
    results,
    factors
) {
  
  lik_df <- data.frame()
  
  for (n_factor in factors) {
    
    factors_str <- glue::glue("{n_factor}_factors")
    
    fastglmpca_28cores_fit <- results[["fastglmpca_28_cores"]][[factors_str]]
    
    fastglmpca_1core_fit <- results[["fastglmpca_1_core"]][[factors_str]]
    
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
      dplyr::mutate(dist_from_best_ll = abs(loglik - best_ll) + 1)
    
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
      scale_y_continuous(
        trans = "log10",
        breaks = 10 ^ seq(0, floor(log10(max(factor_df$dist_from_best_ll)) + 1))
      ) +
      scale_x_continuous(breaks = 0:10) + 
      xlab("Time (hours)") +
      ylab("Dist. from Best Log-Lik") +
      ggtitle(glue::glue("{dataset} {n_factor} Factor Progress")) +
      scale_color_manual(
        values = c(
          "royalblue",
          "darkorange",
          "forestgreen",
          "gold"
        )
      ) +
      theme(panel.border = element_blank(), axis.line = element_line())
    
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
  c(10)
)

droplet_plots_list <- create_plot_list(
  droplet_df,
  "Droplets",
  c(10)
)

pbmc_df <- create_dataset_df(
  pbmc_res_list,
  c(10)
)

pbmc_plots_list <- create_plot_list(
  pbmc_df,
  "PBMC 68k",
  c(10)
)

progress_ll <- ggarrange(
  droplet_plots_list$l10_10, 
  pbmc_plots_list$l10_10, 
  nrow = 1, 
  ncol = 2, 
  common.legend = TRUE, legend = "right",
  labels = c("C", "D")
)

ggarrange(
  pcs, progress_ll, 
  nrow = 2, ncol = 1,
  heights = c(1, .5)
)

