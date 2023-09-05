# Now, I have everything required for the loglik plots.
# It would be great to 

pbmc_res_list <- list()
  
load("~/Downloads/pbmc_68k.RData")

ll_const <- sum(MatrixExtra::mapSparse(counts, lfactorial))

for (model in c("scGBM", "glmpca", "fastglmpca_fit_1_core", "fastglmpca_fit_28_cores")) {
  
  pbmc_res_list[[model]] <- list()
  
  for (factors in c(2, 3, 4, 5, 10, 15, 25)) {
    
    if (model == "scGBM") {
      
      mod_str <- glue::glue(
        "pbmc_scGBM_fit_{factors}_factors_no_beta_infer_10_hrs.rds"
      )
      
    } else if (model == "glmpca") {
      
      mod_str <- glue::glue(
        "pbmc_glmpca_fit_{factors}_factors_avagrad_optimizer_minibatch_stochastic_10_hrs.rds"
      )
      
    } else {
      
      mod_str <- glue::glue(
        "pbmc_{model}_{factors}_factors_10_hrs.rds"
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
      
      pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["loglik"]] <- mod$progress$loglik
      pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["time"]] <- cumsum(mod$progress$time)
      
      if (factors == 2) {
        
        tmp_mod <- fastglmpca:::orthonormalize(mod$U[,3:4], mod$V[,3:4])
        pbmc_res_list[[model]][[glue::glue("{factors}_factors")]][["V"]] <- tmp_mod$V
        
      }
      
    }
    
    
  }
  
}

load("~/Downloads/droplet.RData")
ll_const <- sum(MatrixExtra::mapSparse(counts, lfactorial))

droplets_res_list <- list()

for (model in c("scGBM", "glmpca", "fastglmpca_fit_1_core", "fastglmpca_fit_28_cores")) {
  
  droplets_res_list[[model]] <- list()
  
  for (factors in c(2, 3, 4, 5, 10, 15, 25)) {
    
    if (model == "scGBM") {
      
      mod_str <- glue::glue(
        "droplets_scGBM_fit_{factors}_factors_no_beta_infer_10_hrs.rds"
      )
      
    } else if (model == "glmpca") {
      
      mod_str <- glue::glue(
        "droplets_glmpca_fit_{factors}_factors_avagrad_optimizer_minibatch_stochastic_10_hrs.rds"
      )
      
    } else {
      
      mod_str <- glue::glue(
        "droplets_{model}_{factors}_factors_10_hrs.rds"
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
      
      droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["loglik"]] <- mod$progress$loglik
      droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["time"]] <- cumsum(mod$progress$time)
      
      if (factors == 2) {
        
        tmp_mod <- fastglmpca:::orthonormalize(mod$U[,3:4], mod$V[,3:4])
        droplets_res_list[[model]][[glue::glue("{factors}_factors")]][["V"]] <- tmp_mod$V
        
      }
      
    }
    
    
  }
  
}

save(
  pbmc_res_list, 
  droplets_res_list,
  file = "~/Documents/fastglmpca/inst/analysis/results.RData"
)
