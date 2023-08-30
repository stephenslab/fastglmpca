library(dplyr)
library(ggplot2)

create_dataset_df <- function(
  dataset,
  factors
) {
  
  if (dataset == "droplets") {
    
    load("~/Downloads/droplet.RData")
    
  } else if (dataset == "pbmc") {
    
    load("~/Downloads/pbmc_68k.RData")
    
  }
  
  ll_const <- sum(MatrixExtra::mapSparse(counts, lfactorial))
  
  lik_df <- data.frame()
  
  for (n_factor in factors) {
    
    fastglmpca_28cores_fit <- readr::read_rds(
      glue::glue(
        "{dataset}_fastglmpca_fit_28_cores_{n_factor}_factors_10_hrs.rds"
      )
    )
    
    fastglmpca_1core_fit <- readr::read_rds(
      glue::glue(
        "{dataset}_fastglmpca_fit_1_core_{n_factor}_factors_10_hrs.rds"
      )
    )
    
    glmpca_fit <- readr::read_rds(
      glue::glue(
        "{dataset}_glmpca_fit_{n_factor}_factors_avagrad_optimizer_minibatch_stochastic_10_hrs.rds"
      )
    )
    
    scGBM_fit <- readr::read_rds(
      glue::glue(
        "{dataset}_scGBM_fit_{n_factor}_factors_no_beta_infer_10_hrs.rds"
      )
    )
    
    lik_factor_df <- data.frame(
      factors = rep(
        n_factor, 
        length(fastglmpca_28cores_fit$progress$iter) + length(fastglmpca_1core_fit$progress$iter) +
          length(scGBM_fit$loglik) + length(glmpca_fit$lik)
        ),
      Algorithm = c(
        rep("fastglmpca-28core", length(fastglmpca_28cores_fit$progress$iter)),
        rep("fastglmpca-1core", length(fastglmpca_1core_fit$progress$iter)),
        rep("scGBM", length(scGBM_fit$loglik)),
        rep("glmpca-avagrad-sgd", length(glmpca_fit$lik))
      ),
      loglik = c(
        fastglmpca_28cores_fit$progress$loglik,
        fastglmpca_1core_fit$progress$loglik,
        scGBM_fit$loglik - ll_const,
        glmpca_fit$lik
      ),
      time = c(
        cumsum(fastglmpca_28cores_fit$progress$time),
        cumsum(fastglmpca_1core_fit$progress$time),
        scGBM_fit$time,
        cumsum(c(0, glmpca_fit$time))
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
    
    g_l10 <- ggplot(data = factor_df) +
      geom_point(aes(x = time / 3600, y = dist_from_best_ll, color = Algorithm)) +
      geom_line(aes(x = time / 3600, y = dist_from_best_ll, color = Algorithm)) + 
      scale_y_continuous(trans = "log10") +
      xlab("Time (hours)") +
      ylab("Dist. from Best Log-Lik") + 
      ggtitle(glue::glue("{dataset} {n_factor} Factor")) +
      scale_color_brewer(palette = "Spectral")
    
    #g_std <- ggplot(data = factor_df) +
    #  geom_point(aes(x = time / 3600, y = loglik, color = algo)) +
    #  geom_line(aes(x = time/ 3600, y = loglik, color = algo)) + 
    #  xlab("Time (hours)") +
    #  ylab("Log Likelihood")
    
    l10_list[[glue::glue("l10_{n_factor}")]] <- g_l10
    #std_list[[glue::glue("std_{n_factor}")]] <- g_std
    
  }
  
  return(
    c(
      l10_list
      #std_list
    )
  )
  
}

create_plot_list_std <- function(
    dataset_df,
    dataset,
    factors
) {
  
  #l10_list <- list()
  std_list <- list()
  
  for (n_factor in factors) {
    
    factor_df <- dataset_df %>%
      dplyr::filter(factors == n_factor)
    
    y_upper <- max(factor_df$loglik)
    f_df <- factor_df %>% filter(Algorithm == "fastglmpca-28core")
    y_lower <- min(f_df$loglik)
    
    #g_l10 <- ggplot(data = factor_df) +
    #  geom_point(aes(x = time / 3600, y = dist_from_best_ll, color = Algorithm)) +
    #  geom_line(aes(x = time / 3600, y = dist_from_best_ll, color = Algorithm)) + 
    #  scale_y_continuous(trans = "log10") +
    #  xlab("Time (hours)") +
    #  ylab("Dist. from Best Log-Lik") + 
    #  ggtitle(glue::glue("{dataset} {n_factor} Factor")) +
    #  scale_color_brewer(palette = "Spectral")
    
    g_std <- ggplot(data = factor_df) +
      geom_point(aes(x = time / 3600, y = loglik, color = Algorithm)) +
      geom_line(aes(x = time/ 3600, y = loglik, color = Algorithm)) + 
      xlab("Time (hours)") +
      ylab("Log Likelihood") +
      ylim(y_lower, y_upper) +
      ggtitle(glue::glue("{dataset} {n_factor} Factor")) +
      scale_color_brewer(palette = "Spectral")
    
    #l10_list[[glue::glue("l10_{n_factor}")]] <- g_l10
    std_list[[glue::glue("std_{n_factor}")]] <- g_std
    
  }
  
  return(
    c(
      #l10_list
      std_list
    )
  )
  
}

droplet_df <- create_dataset_df(
  "droplets",
  c(2,3,4,5,10,15,25)
)

droplet_plots_list <- create_plot_list(
    droplet_df,
    "Droplets",
    c(2,3,4,5,10,15,25)
)

droplet_plots_list_std <- create_plot_list_std(
  droplet_df,
  "Droplets",
  c(2,3,4,5,10,15,25)
)

library(ggpubr)

ggarrange(plotlist = droplet_plots_list, nrow = 3, ncol = 3, common.legend = TRUE)

pbmc_df <- create_dataset_df(
  "pbmc",
  c(2,3,4,5,10,15,25)
)

pbmc_plots_list <- create_plot_list(
  pbmc_df,
  "PBMC",
  c(2,3,4,5,10,15,25)
)

library(ggpubr)

ggarrange(plotlist = pbmc_plots_list, nrow = 3, ncol = 3, common.legend = TRUE)
