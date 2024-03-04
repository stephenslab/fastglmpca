load("~/Documents/fastglmpca/inst/analysis/results.RData")

find_first_greater_index <- function(x, y) {
  greater_indices <- which(y > x)
  if (length(greater_indices) > 0) {
    return(min(greater_indices))
  } else {
    return(NA)
  }
}

K_vec <- c(2, 3, 4, 5, 10, 15, 25)

time_until_beats_glmpca_1core <- c()
time_until_beats_glmpca_28core <- c()

time_until_beats_scGBM_1core <- c()
time_until_beats_scGBM_28core <- c()

for (K in K_vec) {
  
  ll_glmpca <- droplets_res_list$glmpca[[
    glue::glue("{K}_factors")
  ]][["loglik"]]
  
  ll_scGBM <- droplets_res_list$scGBM[[
    glue::glue("{K}_factors")
  ]][["loglik"]]
  
  best_ll_glmpca <- tail(ll_glmpca, 1)
  best_ll_scGBM <- tail(ll_scGBM, 1)
  
  ll_fastglmpca_1core <- droplets_res_list$fastglmpca_1_core[[
    glue::glue("{K}_factors")
  ]][["loglik"]]
  
  ll_fastglmpca_28cores <- droplets_res_list$fastglmpca_28_cores[[
    glue::glue("{K}_factors")
  ]][["loglik"]]
  
  idx_1core_glmpca <- find_first_greater_index(
    best_ll_glmpca, 
    ll_fastglmpca_1core
  )
  
  idx_28core_glmpca <- find_first_greater_index(
    best_ll_glmpca, 
    ll_fastglmpca_28cores
  )
  
  idx_1core_scGBM <- find_first_greater_index(
    best_ll_scGBM, 
    ll_fastglmpca_1core
  )
  
  idx_28core_scGBM <- find_first_greater_index(
    best_ll_scGBM, 
    ll_fastglmpca_28cores
  )
  
  if(is.na(idx_1core_glmpca)) {
    
    time_until_beats_glmpca_1core <- c(
      time_until_beats_glmpca_1core, 
      Inf
    )
    
  } else {
    
    time_until_beats_glmpca_1core <- c(
      time_until_beats_glmpca_1core, 
      droplets_res_list$fastglmpca_1_core[[
        glue::glue("{K}_factors")
      ]][["time"]][idx_1core_glmpca] / 3600
    )
    
  }
  
  if(is.na(idx_28core_glmpca)) {
    
    time_until_beats_glmpca_28core <- c(
      time_until_beats_glmpca_28core, 
      Inf
    )
    
  } else {
    
    time_until_beats_glmpca_28core <- c(
      time_until_beats_glmpca_28core, 
      droplets_res_list$fastglmpca_28_cores[[
        glue::glue("{K}_factors")
      ]][["time"]][idx_28core_glmpca] / 3600
    )
    
  }
  
  if(is.na(idx_1core_scGBM)) {
    
    time_until_beats_scGBM_1core <- c(
      time_until_beats_scGBM_1core, 
      Inf
    )
    
  } else {
    
    time_until_beats_scGBM_1core <- c(
      time_until_beats_scGBM_1core, 
      droplets_res_list$fastglmpca_1_core[[
        glue::glue("{K}_factors")
      ]][["time"]][idx_1core_scGBM] / 3600
    )
    
  }
  
  if(is.na(idx_28core_scGBM)) {
    
    time_until_beats_scGBM_28core <- c(
      time_until_beats_scGBM_28core, 
      Inf
    )
    
  } else {
    
    time_until_beats_scGBM_28core <- c(
      time_until_beats_scGBM_28core, 
      droplets_res_list$fastglmpca_28_cores[[
        glue::glue("{K}_factors")
      ]][["time"]][idx_28core_scGBM] / 3600
    )
    
  }
  
}

df_droplets <- data.frame(
  time_until_beats_glmpca_1core = time_until_beats_glmpca_1core,
  time_until_beats_glmpca_28core = time_until_beats_glmpca_28core,
  time_until_beats_scGBM_1core = time_until_beats_scGBM_1core,
  time_until_beats_scGBM_28core = time_until_beats_scGBM_28core
)

library(dplyr)

df_droplets <- df_droplets %>%
  dplyr::mutate(
    time_until_beats_glmpca_1core = if_else(
      is.infinite(time_until_beats_glmpca_1core),
      "> 10",
      as.character(round(time_until_beats_glmpca_1core, 1))
    ),
    time_until_beats_glmpca_28core = if_else(
      is.infinite(time_until_beats_glmpca_28core),
      "> 10",
      as.character(round(time_until_beats_glmpca_28core, 2))
    ),
    time_until_beats_scGBM_1core = as.character(round(time_until_beats_scGBM_1core, 2)),
    time_until_beats_scGBM_28core = as.character(round(time_until_beats_scGBM_28core, 3))
  )

time_until_beats_glmpca_1core <- c()
time_until_beats_glmpca_28core <- c()

time_until_beats_scGBM_1core <- c()
time_until_beats_scGBM_28core <- c()

for (K in K_vec) {
  
  ll_glmpca <- pbmc_res_list$glmpca[[
    glue::glue("{K}_factors")
  ]][["loglik"]]
  
  ll_scGBM <- pbmc_res_list$scGBM[[
    glue::glue("{K}_factors")
  ]][["loglik"]]
  
  best_ll_glmpca <- tail(ll_glmpca, 1)
  best_ll_scGBM <- tail(ll_scGBM, 1)
  
  ll_fastglmpca_1core <- pbmc_res_list$fastglmpca_1_core[[
    glue::glue("{K}_factors")
  ]][["loglik"]]
  
  ll_fastglmpca_28cores <- pbmc_res_list$fastglmpca_28_cores[[
    glue::glue("{K}_factors")
  ]][["loglik"]]
  
  idx_1core_glmpca <- find_first_greater_index(
    best_ll_glmpca, 
    ll_fastglmpca_1core
  )
  
  idx_28core_glmpca <- find_first_greater_index(
    best_ll_glmpca, 
    ll_fastglmpca_28cores
  )
  
  idx_1core_scGBM <- find_first_greater_index(
    best_ll_scGBM, 
    ll_fastglmpca_1core
  )
  
  idx_28core_scGBM <- find_first_greater_index(
    best_ll_scGBM, 
    ll_fastglmpca_28cores
  )
  
  if(is.na(idx_1core_glmpca)) {
    
    time_until_beats_glmpca_1core <- c(
      time_until_beats_glmpca_1core, 
      Inf
    )
    
  } else {
    
    time_until_beats_glmpca_1core <- c(
      time_until_beats_glmpca_1core, 
      pbmc_res_list$fastglmpca_1_core[[
        glue::glue("{K}_factors")
      ]][["time"]][idx_1core_glmpca] / 3600
    )
    
  }
  
  if(is.na(idx_28core_glmpca)) {
    
    time_until_beats_glmpca_28core <- c(
      time_until_beats_glmpca_28core, 
      Inf
    )
    
  } else {
    
    time_until_beats_glmpca_28core <- c(
      time_until_beats_glmpca_28core, 
      pbmc_res_list$fastglmpca_28_cores[[
        glue::glue("{K}_factors")
      ]][["time"]][idx_28core_glmpca] / 3600
    )
    
  }
  
  if(is.na(idx_1core_scGBM)) {
    
    time_until_beats_scGBM_1core <- c(
      time_until_beats_scGBM_1core, 
      Inf
    )
    
  } else {
    
    time_until_beats_scGBM_1core <- c(
      time_until_beats_scGBM_1core, 
      pbmc_res_list$fastglmpca_1_core[[
        glue::glue("{K}_factors")
      ]][["time"]][idx_1core_scGBM] / 3600
    )
    
  }
  
  if(is.na(idx_28core_scGBM)) {
    
    time_until_beats_scGBM_28core <- c(
      time_until_beats_scGBM_28core, 
      Inf
    )
    
  } else {
    
    time_until_beats_scGBM_28core <- c(
      time_until_beats_scGBM_28core, 
      pbmc_res_list$fastglmpca_28_cores[[
        glue::glue("{K}_factors")
      ]][["time"]][idx_28core_scGBM] / 3600
    )
    
  }
  
}

df_pbmc <- data.frame(
  time_until_beats_glmpca_1core = time_until_beats_glmpca_1core,
  time_until_beats_glmpca_28core = time_until_beats_glmpca_28core,
  time_until_beats_scGBM_1core = time_until_beats_scGBM_1core,
  time_until_beats_scGBM_28core = time_until_beats_scGBM_28core
)

library(dplyr)

df_pbmc <- df_pbmc %>%
  dplyr::mutate(
    time_until_beats_glmpca_1core = if_else(
      is.infinite(time_until_beats_glmpca_1core),
      "> 10",
      as.character(round(time_until_beats_glmpca_1core, 1))
    ),
    time_until_beats_glmpca_28core = if_else(
      is.infinite(time_until_beats_glmpca_28core),
      "> 10",
      as.character(round(time_until_beats_glmpca_28core, 2))
    ),
    time_until_beats_scGBM_1core = as.character(round(time_until_beats_scGBM_1core, 2)),
    time_until_beats_scGBM_28core = as.character(round(time_until_beats_scGBM_28core, 3))
  )

rownames(df_pbmc) <- c("K = 2", "K = 3", "K = 4", "K = 5", "K = 10", "K = 15", "K = 25")
rownames(df_droplets) <- c("K = 2", "K = 3", "K = 4", "K = 5", "K = 10", "K = 15", "K = 25")


