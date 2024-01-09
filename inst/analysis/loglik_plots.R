library(dplyr)
library(ggplot2)

load("~/Documents/updated_fastglmpca/fastglmpca/inst/analysis/results.RData")

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
      dplyr::mutate(dist_from_best_ll = abs(loglik - best_ll + .1))

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
      ylab("Dist. from Best Log-Lik") +
      ggtitle(glue::glue("{dataset} {n_factor} Factor")) +
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
  c(2,3,4,5,10,15,25)
)

droplet_plots_list <- create_plot_list(
    droplet_df,
    "Droplets",
    c(2,3,4,5,10,15,25)
)

library(ggpubr)

ggarrange(
  plotlist = droplet_plots_list,
  nrow = 3,
  ncol = 3,
  common.legend = TRUE,
  legend = "right"
  )

pbmc_df <- create_dataset_df(
  pbmc_res_list,
  c(2,3,4,5,10,15,25)
)

pbmc_plots_list <- create_plot_list(
  pbmc_df,
  "PBMC",
  c(2,3,4,5,10,15,25)
)

library(ggpubr)

ggarrange(
  plotlist = pbmc_plots_list,
  nrow = 3,
  ncol = 3,
  common.legend = TRUE,
  legend = "right"
  )
