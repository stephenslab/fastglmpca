
pbmc_res_list <- list()

load("~/Downloads/pbmc_68k.RData")

ll_const <- sum(MatrixExtra::mapSparse(counts, lfactorial))

iter_run_28core <- c(2200, 1500, 1050, 850, 450, 150, 200)
iter_run_1core <- c(85, 55, 50, 40, 25, 18, 12)

factor_vec <- c(2, 3, 4, 5, 10, 15, 25)

for (model in c("scGBM", "glmpca", "fastglmpca_1_core", "fastglmpca_28_cores")) {

  pbmc_res_list[[model]] <- list()

  for (i in 1:length(factor_vec)) {

    factors <- factor_vec[i]

    if (model == "scGBM") {

      mod_str <- glue::glue(
        "pbmc_scGBM_fit_{factors}_factors_no_beta_infer_10_hrs.rds"
      )

    } else if (model == "glmpca") {

      mod_str <- glue::glue(
        "pbmc_glmpca_fit_{factors}_factors_10_hrs_avagrad_optimizer_minibatch_stochastic_dec_23.rds"
      )

    } else if (model == "fastglmpca_1_core") {

      mod_str <- glue::glue(
        "pbmc_fastglmpca_fit_{factors}_factors_{iter_run_1core[i]}_iter_1_core_dec_23.rds"
      )

    } else if (model == "fastglmpca_28_cores") {

      mod_str <- glue::glue(
        "pbmc_fastglmpca_fit_{factors}_factors_{iter_run_28core[i]}_iter_28_cores_dec_23.rds"
      )

    }

    mod <- readr::read_rds(mod_str)

    if (model == "scGBM") {

      pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["loglik"]] <- mod$loglik - ll_const
      pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["time"]] <- mod$time - mod$time[1]

      if (factors == 2) {

        pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["V"]] <- mod$V

      }

    } else if (model == "glmpca") {

      pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["loglik"]] <- mod$lik
      pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["time"]] <- cumsum(c(0, mod$time))

      if (factors == 2) {

        tmp_mod <- fastglmpca:::orthonormalize(
          as.matrix(mod$loadings), as.matrix(mod$factors)
        )

        pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["V"]] <- tmp_mod$V

      }

    } else {

      if (max(cumsum(mod$progress$time)) > 60 * 60 * 10) {

        idx_10hr <- min(which(cumsum(mod$progress$time) > (60 * 60 * 10))) - 1

      } else {

        idx_10hr <- length(mod$progress$time)

      }

      pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["loglik"]] <- mod$progress$loglik[1:idx_10hr]
      pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["time"]] <- cumsum(mod$progress$time)[1:idx_10hr]

      if (factors == 2 && model == "fastglmpca_28_cores") {

        pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["V"]] <- mod$V

      }

    }

  }

}


load("~/Downloads/droplet.RData")
ll_const <- sum(MatrixExtra::mapSparse(counts, lfactorial))

iter_run_28core <- c(28005, 23005, 14005, 12505, 5105, 3405, 2105)
iter_run_1core <- c(1400, 700, 550, 400, 225, 150, 105)

factor_vec <- c(2, 3, 4, 5, 10, 15, 25)

droplets_res_list <- list()

for (model in c("scGBM", "glmpca", "fastglmpca_1_core", "fastglmpca_28_cores")) {

  droplets_res_list[[model]] <- list()

  for (i in 1:length(factor_vec)) {

    factors <- factor_vec[i]

    if (model == "scGBM") {

      mod_str <- glue::glue(
        "droplets_scGBM_fit_{factors}_factors_no_beta_infer_10_hrs.rds"
      )

    } else if (model == "glmpca") {

      mod_str <- glue::glue(
        "droplets_glmpca_fit_{factors}_factors_10_hrs_avagrad_optimizer_minibatch_stochastic_dec_23.rds"
      )

    } else if (model == "fastglmpca_1_core") {

      mod_str <- glue::glue(
        "droplets_fastglmpca_fit_{factors}_factors_{iter_run_1core[i]}_iter_1_core_dec_23.rds"
      )

    } else if (model == "fastglmpca_28_cores") {

      mod_str <- glue::glue(
        "droplets_fastglmpca_fit_{factors}_factors_{iter_run_28core[i]}_iter_28_cores_dec_23.rds"
      )

    }

    mod <- readr::read_rds(mod_str)

    if (model == "scGBM") {

      droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["loglik"]] <- mod$loglik - ll_const
      droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["time"]] <- mod$time - mod$time[1]

      if (factors == 2) {

        droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["V"]] <- mod$V

      }

    } else if (model == "glmpca") {

      droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["loglik"]] <- mod$lik
      droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["time"]] <- cumsum(c(0, mod$time))

      if (factors == 2) {

        tmp_mod <- fastglmpca:::orthonormalize(
          as.matrix(mod$loadings), as.matrix(mod$factors)
        )

        droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["V"]] <- tmp_mod$V

      }

    } else {

      # First, need to make sure that I'm only taking 10 hours worth of data
      if (max(cumsum(mod$progress$time)) > 60 * 60 * 10) {

        idx_10hr <- min(which(cumsum(mod$progress$time) > (60 * 60 * 10))) - 1

      } else {

        idx_10hr <- length(mod$progress$time)

      }

      droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["loglik"]] <- mod$progress$loglik[1:idx_10hr]
      droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["time"]] <- cumsum(mod$progress$time)[1:idx_10hr]

      if (factors == 2 && model == "fastglmpca_28_cores") {

        droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["V"]] <- mod$V

      }

    }

  }

}

save(
  pbmc_res_list,
  droplets_res_list,
  file = "~/Documents/updated_fastglmpca/fastglmpca/inst/analysis/results.RData"
)
